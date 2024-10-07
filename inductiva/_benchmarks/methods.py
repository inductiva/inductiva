"""
Methods to interact with Benchmarks.
"""
from collections import defaultdict
import datetime
import logging
import json

from typing import Any, Callable, Dict, List, Type, Union
from tqdm import tqdm

import inductiva
from inductiva.client import models
from inductiva.resources.machines_base import BaseMachineGroup
from inductiva.templating.manager import TemplateManager
from inductiva.resources.machines import MachineGroup
from inductiva.simulators.simulator import Simulator
from inductiva.projects.project import Project
from inductiva.resources import machine_groups
from inductiva.tasks.task import Task

from ..localization import translator as __


def _analyze_task(
    task: Task,
    report: Dict[str, Any],
    skipped_tasks: List[Task],
):
    """Analyzes a task.

    This method analyzes a task and adds the information to the report. If the
    task is not successful, it is added to the skipped_tasks list.

    Args:
        task (Task): Task to analyze.
        report (Dict[str, Any]): Report to add the information to.
        skipped_tasks (List[Task]): List of tasks that were skipped.
    Returns:
        This method does not return anything. It modifies some parameters
        instead of returning multiple values.
    """
    task_status = task.get_status()
    if task_status != models.TaskStatusCode.SUCCESS:
        skipped_tasks.append(task)
        return
    task_info = task.info
    a, b, *_ = task_info.executer.vm_name.split("-")
    mg_name = a + "-" + b
    resource = machine_groups.get_by_name(mg_name)
    vm_type = task_info.executer.vm_type

    if vm_type not in report["resources"]:
        report["resources"][vm_type] = {
            "cost_per_hour": resource.estimate_cloud_cost(verbose=False),
            "type": resource.__class__.__name__,
            "provider": task_info.executer.host_type,
            "machine_count": resource.num_machines,
            "tasks": []
        }
    metadata_filename = "input.json"
    inputs_path = task.download_inputs(filenames=[metadata_filename])
    metadata_path = inputs_path.joinpath(metadata_filename)
    with open(metadata_path, mode="r", encoding="utf-8") as metadata_file:
        metadata_content = json.load(metadata_file)
    report["resources"][vm_type]["tasks"].append({
        "task_id": task_info.task_id,
        "info": task_info.to_dict(),
        "metadata": metadata_content,
    })


def analyze(name: str, output_path: str = None) -> Dict:
    """Analyzes a benchmark.
    This method creates a dictionary with all the relevant information about
    the benchmark.
    Args:
        name (str): Name of the benchmark.
        output_path (str): Path to the output file.
            If None it will not create an output file.
    Returns:
        Dict: Dictionary with the information about the benchmark.
    """
    logging.info("Starting analysis of benchmark %s", name)
    project = Project(name)

    tasks = project.get_tasks()

    report = {"name": name, "resources": {}}

    skipped_tasks = []
    for task in tqdm(tasks, desc="Analyzing tasks"):
        _analyze_task(task, report, skipped_tasks)

    for task in skipped_tasks:
        logging.info("Task %s was skipped due to status %s", task.info.task_id,
                     task.info.status)
    if output_path:
        logging.info("Writing report to %s/%s.json", output_path, name)
        with open(f"{output_path}/{name}.json", "w",
                  encoding="utf-8") as output_file:
            output_file.write(json.dumps(report, indent=4))

    return report


def _tasks_by_vm_type(project: Project) -> Dict[str, List[Task]]:
    """Returns a dictionary with the resources and the tasks ran on them.
    Args:
        project (Project): Project to get the resources and tasks from.
    Returns:
        Dict[str, List[Task]]: Dictionary with the resources as keys and the
            tasks ran on them as values.
        Example:.
        { "n1-standard-4": [task1, task2], "n1-standard-8": [task3]}
    """
    tasks = project.get_tasks()
    resources = defaultdict(list)
    for task in tasks:
        task_info = task.get_info()
        if task_info.executer is None:
            raise RuntimeError(__("tasks_not_executed_yet"))

        resource = task_info.executer.vm_type
        resources[resource].append(task)
    return resources


def _compute_tasks_to_run(machines_to_run: List[str],
                          already_ran: Dict[str, List[Task]],
                          intended_replicas: int) -> Dict[str, int]:
    """Computes the number of tasks that need to be run for each resource.
    Args:
        machines_to_run (List[str]): List of machines to run the simulation on.
        already_ran (Dict[str, List[Task]]): Dictionary with the resources as
            keys and the tasks ran on them as values.
        intended_replicas (int): Number of replicas to run on each resource.
    Returns:
        Dict[str, int]: Dictionary with the resources as keys and the number of
            tasks that need to be run on them as values.
        Example: { "n1-standard-4": 2, "n1-standard-8": 1}
    """
    already_ran = defaultdict(list, already_ran)
    return {m: intended_replicas - len(already_ran[m]) for m in machines_to_run}


def _print_info(already_ran: Dict[str, List[Task]], tasks_to_run: Dict[str,
                                                                       int]):
    """
    Prints the information about the tasks that have already been run and the
    tasks that need to be run.
    """

    print()
    print(__("tasks_ran"))
    for key, value in already_ran.items():
        print(f"--->{key}: {len(value)}")
    print()
    print(__("tasks_to_run"))
    for key, value in tasks_to_run.items():
        print(f"--->{key}: {value}")
    print()


def _run_replica(simulator_args: Dict[str, Union[Callable, Any]],\
                 input_args: Dict[str, Union[Callable,Any]],\
                 resource: BaseMachineGroup,\
                 simulator: Simulator,\
                 input_files: str,\
                 replica: int):
    """Runs a replica of the benchmark.

    Args:
        simulator_args (Dict[str, Union[Callable, Any]]): Arguments for the
            simulator.
        input_args (Dict[str, Union[Callable, Any]]): Arguments for the
            input files.
        resource (BaseMachineGroup): Resource to run the benchmark on.
        simulator (Simulator): Simulator to run.
        input_files (str): Path to the input files.
        replica (int): Number of the replica.
    """
    # Each replica can run with different arguments
    simulator_args_current = _render_dict(simulator_args, resource)
    input_args_current = _render_dict(input_args, resource)
    print(f"Simulation arguments {simulator_args_current}")
    # input_files either come from the template manager or
    # are passed directly.
    simulator_args_current["input_dir"] = input_files
    simulator_args_current["on"] = resource

    if input_args_current:
        target_dir = f"benchmark_inputs/{resource.machine_type}_{replica}"
        TemplateManager.render_dir(source_dir=input_files,
                                   target_dir=target_dir,
                                   overwrite=False,
                                   **input_args)
        simulator_args_current["input_dir"] = target_dir

    _ = simulator.run(**simulator_args_current)


def _run(project: Project,
         input_files: str,
         replicas: int,
         machines: List[str],
         simulator: Simulator,
         input_args: Dict[str, Union[Callable, Any]] = None,
         machine_args: Dict[str, Union[Callable, Any]] = None,
         machine_class: Type[BaseMachineGroup] = inductiva.resources.
         MachineGroup,
         simulator_args: Dict[str, Union[Callable, Any]] = None):
    """Runs a benchmark.
    This method runs a benchmark using the provided parameters.
    Args:
        project (Project): Project to run the benchmark on.
    """
    # create dict with resource type: list of tasks ran on that resource
    # Example
    # { "n1-standard-4": [task1, task2], "n1-standard-8": [task3]}
    current_project_tasks = _tasks_by_vm_type(project)
    to_run = _compute_tasks_to_run(machines, current_project_tasks, replicas)

    print()
    if not to_run:
        print(__("all_tasks_run"))
        project.close()
        return

    _print_info(current_project_tasks, to_run)

    for machine, replicas_to_run in to_run.items():

        machine_args_current = _render_dict(machine_args, machine)

        print(f"\nCreating machine group {machine} with {machine_args_current}")

        resource = machine_class(machine_type=machine,
                                 max_idle_time=datetime.timedelta(minutes=1),
                                 **machine_args_current)

        print(f"\nRunning {simulator}")
        if resource is not None:
            resource.start(wait_on_pending_quota=True)

        for replica in range(replicas_to_run, 0, -1):
            _run_replica(simulator_args, input_args, resource, simulator,
                         input_files, replica)


def run(name: str,
        input_files: str,
        replicas: int,
        machines: List[str],
        simulator: Simulator,
        append: bool = False,
        input_args: Dict[str, Union[Callable, Any]] = None,
        machine_class: Type[BaseMachineGroup] = MachineGroup,
        machine_args: Dict[str, Union[Callable, Any]] = None,
        simulator_args: Dict[str, Union[Callable, Any]] = None):
    """
    Runs a benchmark.

    This method creates a project and runs a benchmark using the provided
    parameters. It checks if a benchmark with the given name already exists. If
    it does, the user is prompted to either append tasks to the existing
    benchmark or cancel the operation. The method waits for quotas to be
    available before starting the necessary resources. If the benchmark is
    already completed, no further action is taken. If the benchmark is
    incomplete, only the remaining tasks are executed.

    input_args, machine_args and simulator_args can be dictionaries where the
    values can be callables that will be called with the resource object as
    parameter (or str machine_type when it comes to machine_args). This is
    useful when you want to calculate the values based on the machine_group or
    machine_type.
    Example:
        def data_disk_gb(name):
            #makes the data disk size relative to the number of vcpus
            vcpu = int(name.split("-")[2])
            return vcpu+10
        machine_args={"data_disk_gb": data_disk_gb}

    Args:
            name (str): Name of the benchmark.
            input_files (str): Path to the input files.
            replicas (int): Number of replicas to run on each resource.
            machines (List[str]): List of machines to run the simulation on.
                The machine in this case is just a string with the machine type.
                Ex: ["c2-standard-4", "c2d-standard-4"]
            simulator (Simulator): Simulator to run.
            append (bool): If True, appends the tasks to a benchmark regardless
            if that benchmark is new or not. If False, will prompt the user if
            the benchmark already exists.
            machine_class (BaseMachineGroup): Class to use to create the machine
                group. If not given, MachineGroup will be used.
            input_args (Dict[str, Union[Callable, Any]]): Arguments for the
                input files.
            machine_args (Dict[str, Union[Callable, Any]]): Arguments for the
                machines.
            simulator_args (Dict[str, Union[Callable, Any]]): Arguments for the
                simulator.
    """

    project = Project(name, append=True)
    if project.num_tasks > 0 and not append:
        raise RuntimeError(
            __("benchmark-already-exists", name, project.num_tasks))

    with project:
        _run(project, input_files, replicas, machines, simulator, input_args,
             machine_args, machine_class, simulator_args)


def _render_dict(dictionary: Dict[str, Union[Callable, Any]], argument: Any):
    """Replaces callables from a dictionary.
    This function replaces callables with the result of calling the callable
    with the argument passed.
    Args:
        dictionary (Dict[str, Union[Callable, Any]]): Dictionary to replace
            callables from.
        argument: Argument to pass to the callables.
    Returns:
        Dict[str, Any]: Dictionary with callables replaced.
    """
    if dictionary:
        return {
            k: v(argument) if callable(v) else v for k, v in dictionary.items()
        }
    return {}
