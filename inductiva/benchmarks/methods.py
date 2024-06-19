""" 
Methods to interact with Benchmarks.
"""
from collections import defaultdict
import datetime
import time
import json

from typing import Any, Callable, Dict, List, Union

from inductiva.resources.machines_base import BaseMachineGroup
from inductiva.templating.manager import TemplateManager
from inductiva.resources.machines import MachineGroup
from inductiva.simulators.simulator import Simulator
from inductiva.projects.project import Project
from inductiva.resources import machine_groups
from inductiva.tasks.task import Task
from inductiva import types, users

from ..localization import translator as __


def analyze(name: str,
            metadata_path: types.Path,
            ouput_path: types.Path = None) -> Dict:
    """Analyzes a benchmark.
    This method creates a dictionary with all the relevant information about
    the benchmark.
    Args:
        name (str): Name of the benchmark.
        metadata_path (types.Path): Path to the metadata file.
        ouput_path (types.Path): Path to the output file.
            If it will not create an output file.
    Returns:
        Dict: Dictionary with the information about the benchmark.
    """

    project = Project(name, append=True)
    project.open()

    tasks = project.get_tasks()

    report = {"name": name, "resources": {}}

    #reads the metadata file and creates a dictionary with the task_id as key
    #and the metadata as value
    metadata_tasks = {}
    with open(metadata_path, "r", encoding="utf-8") as metadata_file:
        lines = metadata_file.readlines()
        for line in lines:
            temp_task = json.loads(line)
            metadata_tasks[temp_task["task_id"]] = temp_task

    for task in tasks:

        if task.get_status() != "success":
            continue
        task_info = task.info
        mg_name = task_info.executer.vm_name.split(
            "-")[0] + "-" + task_info.executer.vm_name.split("-")[1]
        execution_time = task.get_computation_time(cached=True)
        resource = machine_groups.get_by_name(mg_name)
        vm_type = task_info.executer.vm_type

        if vm_type not in report["resources"]:
            report["resources"][vm_type] = {
                "cost_per_hour":
                    resource.estimate_cloud_cost(verbose=False),
                "type":
                    type(resource).__name__,
                "provider":
                    task_info.executer.host_type,
                "machine_count":
                    resource.num_machines,
                "tasks": [{
                    "task_id": task_info.task_id,
                    "start_time": task_info.start_time,
                    "submit_time": task_info.create_time,
                    "execution_time_seconds": execution_time,
                    "metadata": metadata_tasks[task_info.task_id]
                }]
            }
        else:
            report["resources"][vm_type]["tasks"].append({
                "task_id": task_info.task_id,
                "start_time": task_info.start_time,
                "submit_time": task_info.create_time,
                "execution_time_seconds": execution_time,
                "metadata": metadata_tasks[task_info.task_id]
            })
    if ouput_path:
        with open(f"{ouput_path}/{name}.json", "w",
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


def _compute_tasks_to_run(machines_to_run: List[Dict[str, Any]],
                          already_ran: Dict[str, List[Task]],
                          intended_replicas: int) -> Dict[str, int]:
    """Computes the number of tasks that need to be run for each resource.
    Args:
        machines_to_run (List[Dict[str, Any]]): List of machines to run the
            simulation on.
        already_ran (Dict[str, List[Task]]): Dictionary with the resources as
            keys and the tasks ran on them as values.
        intended_replicas (int): Number of replicas to run on each resource.
    Returns:
        Dict[str, int]: Dictionary with the resources as keys and the number of
            tasks that need to be run on them as values.
        Example: { "n1-standard-4": 2, "n1-standard-8": 1}
    """
    to_run = {}
    for machine in machines_to_run:

        # if the machine has not been ran yet, run the intended replicas
        if machine["machine_type"] not in already_ran:
            to_run[machine["machine_type"]] = intended_replicas
        # if the machine has been ran, run the remaining replicas
        elif len(already_ran[machine["machine_type"]]) < intended_replicas:
            to_run[machine["machine_type"]] = intended_replicas - len(
                already_ran[machine["machine_type"]])
    return to_run


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


def run(name: str,
        input_files: types.Path,
        replicas: int,
        machines: List[Dict[str, Any]],
        simulator: Simulator,
        append: bool = False,
        machine_class: BaseMachineGroup = None,
        input_args: Dict[str, Union[Callable, Any]] = None,
        machine_args: Dict[str, Union[Callable, Any]] = None,
        simulator_args: Dict[str, Union[Callable, Any]] = None):
    """Runs a benchmark.
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
            input_files (types.Path): Path to the input files.
            replicas (int): Number of replicas to run on each resource.
            machines (List[Dict[str, Any]]): List of machines to run the
                simulation on.
                The machine in this case is a dict with a lot of information.
                It is mandatory to have the machine_type key and value.
                Usually used in conjunction with inductiva.resources.query
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
    project.open()

    tasks = project.get_tasks()

    if len(tasks) > 0 and not append:
        print(__("benchmark-already-exists", name))
        response = input()
        if response in ("no", "n"):
            print(__("operation_canceled"))
            project.close()
            return
        elif response in ("yes", "ye", "y"):
            print(__("appending_tasks"))

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

    if input_args is not None:
        template_manager = TemplateManager(template_dir=input_files)
        base_path = template_manager.get_root_dir()
    for machine, current_replicas in to_run.items():

        machine_args_current = _replace_callable_from_dict(
            machine_args, machine)

        print(f"\nCreating machine group {machine} with {machine_args_current}")

        # Use MachineGroup as default value for machine_class
        machine_class = machine_class or MachineGroup

        resource = machine_class(machine_type=machine,
                                 max_idle_time=datetime.timedelta(minutes=1),
                                 **machine_args_current)

        simulator_args_current = _replace_callable_from_dict(
            simulator_args, resource)
        input_args_current = _replace_callable_from_dict(input_args, resource)

        if input_args_current:

            template_manager.set_root_dir(f"{base_path}/{name}")
            input_files = template_manager.render_dir(**input_args_current)

        # input_files either come from the template manager or
        # are passed directly
        simulator_args_current["input_dir"] = input_files
        simulator_args_current["on"] = resource

        print(f"\nRunning {simulator}"
              f"with {simulator_args_current}")
        if resource is not None:
            while not _can_start_resource(resource):
                print("This machine will exceed the current quotas.\n"
                      "Going to sleep and trying again later.")
                time.sleep(60)

            resource.start()

        for _ in range(current_replicas):
            _ = simulator.run(**simulator_args_current)

    project.close()


def _replace_callable_from_dict(dictionary: Dict[str, Union[Callable, Any]],
                                argument: Any):
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
    new_dict = {}
    if dictionary is not None:
        for key, value in dictionary.items():
            if callable(value):
                new_dict[key] = value(argument)
            else:
                new_dict[key] = value
    return new_dict


def _can_start_resource(resource: BaseMachineGroup) -> bool:
    """Check if the resource can be started.
        
        This method checks if the resource can be started by checking
        the available quotas and resource usage.

        returns:
            bool: True if the resource can be started, False otherwise.
        """
    quotas = users.get_quotas()
    cost_in_use = quotas["cost_per_hour"]["in_use"]
    vcpu_in_use = quotas["total_num_vcpus"]["in_use"]
    machines_in_use = quotas["total_num_machines"]["in_use"]

    cost_max = quotas["cost_per_hour"]["max_allowed"]
    vcpu_max = quotas["total_num_vcpus"]["max_allowed"]
    machines_max = quotas["total_num_machines"]["max_allowed"]

    current_vcpu = int(resource.machine_type.split("-")[2])

    if (cost_in_use + resource.estimate_cloud_cost(verbose=False) > cost_max or
            vcpu_in_use + current_vcpu > vcpu_max or
            machines_in_use + 1 > machines_max):
        return False
    return True
