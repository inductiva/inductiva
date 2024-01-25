"""Tasks CLI subcommands."""
from inductiva import tasks, utils
from inductiva import _cli
from inductiva.utils import format_utils


def get_task_generic_info(list_of_tasks: list) -> tuple[list, list]:
    columns = [
        "ID", "Simulator", "Status", "Submitted", "Started", "Computation Time",
        "Resource Type"
    ]
    rows = []

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
            if executer["n_mpi_hosts"] > 1:
                resource_type += f" x{executer['n_mpi_hosts']}"

        row = [
            task.id,
            task.get_simulator_name(),
            task.get_status(),
            info.get("input_submit_time", None),
            info.get("start_time", None),
            execution_time,
            resource_type,
        ]
        rows.append(row)

    return rows, columns


def print_tasks_generic_info(rows: list, columns: list):
    formatters = {
        "Submitted": format_utils.datetime_formatter,
        "Started": format_utils.datetime_formatter
    }

    print(
        utils.format_utils.get_tabular_str(
            rows,
            columns,
            formatters=formatters,
        ))


def list_tasks(args):
    """List tasks."""

    list_of_tasks = tasks.list(last_n=args.last_n)

    rows, columns = get_task_generic_info(list_of_tasks)

    print_tasks_generic_info(rows, columns)


def list_task_by_id(args):
    """List a task with a specific ID."""
    t = tasks.Task(args.task_id)

    rows, columns = get_task_generic_info([t])

    print_tasks_generic_info(rows, columns)


def list_command(args):

    if args.task_id is not None:
        list_task_by_id(args)
    else:
        number_of_tasks = 5 if args.last_n is None else args.last_n
        args.last_n = number_of_tasks
        list_tasks(args)


def register_tasks_cli(parser):
    _cli.utils.show_help_msg(parser)

    subparsers = parser.add_subparsers()

    list_subparser = subparsers.add_parser("list", help="List tasks")

    group = list_subparser.add_mutually_exclusive_group()

    group.add_argument("-n",
                       "--last-n",
                       type=int,
                       help="List last N tasks. Default: 5.")

    group.add_argument("-tid",
                       "--task-id",
                       type=str,
                       help="List a task with a specific ID.")

    list_subparser.set_defaults(func=list_command)
