"""Tasks CLI subcommands."""
from inductiva import tasks, utils
from inductiva import _cli
from inductiva.utils import format_utils
from typing import Tuple
from collections import defaultdict


def get_task_generic_info(list_of_tasks: list) -> Mapping[str, List[Any]]:

    table = defaultdict(list)

    for task in list_of_tasks:
        info = task.get_info()
        status = task.get_status()

        computation_end_time = info.get("computation_end_time", None)

        execution_time = task.get_computation_time(fail_if_running=False)
        if execution_time is not None:
            if computation_end_time is None:
                if status in ["started", "submitted"]:
                    execution_time = f"*{execution_time}"
                else:
                    execution_time = "n/a"

        executer = info["executer"]
        if executer is None:
            resource_type = None
        else:
            resource_type = executer["vm_type"]
            n_mpi_hosts = executer["n_mpi_hosts"]
            if n_mpi_hosts > 1:
                resource_type += f" x{n_mpi_hosts}"
        table["ID"].append(task.id)
        table["Simulator"].append(task.get_simulator_name())
        table["Status"].append(status)
        table["Submitted"].append(info.get("input_submit_time", None))
        table["Started"].append(info.get("start_time", None))
        table["Computation Time"].append(execution_time)
        table["Resource Type"].append(resource_type)

    return table


def list_tasks(args):
    """
    List tasks based on the flags used.
    It can print the last N tasks or a single task with a specific ID
    """
    if args.task_id is not None:
        task_list = [tasks.Task(args.task_id)]

    else:
        last_n = 5 if args.last_n is None else args.last_n

        task_list = tasks.list(last_n=last_n)

    table = get_task_generic_info(task_list)

    formatters = {
        "Submitted": format_utils.datetime_formatter,
        "Started": format_utils.datetime_formatter
    }

    print(utils.format_utils.get_tabular_str(
        table,
        formatters=formatters,
    ))


def register_tasks_cli(parser):
    _cli.utils.show_help_msg(parser)

    subparsers = parser.add_subparsers()

    list_subparser = subparsers.add_parser("list", help="List tasks")

    group = list_subparser.add_mutually_exclusive_group()

    group.add_argument("-n",
                       "--last-n",
                       type=int,
                       help="List last N tasks. Default: 5.")

    group.add_argument("-id",
                       "--task-id",
                       type=str,
                       help="List a task with a specific ID.")

    list_subparser.set_defaults(func=list_command)
