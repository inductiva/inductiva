"""Tasks CLI subcommands."""
from inductiva import tasks, utils
from inductiva import _cli
from inductiva.utils import format_utils
from typing import Sequence


def _data_of_tasks(tasks: Sequence["inductiva.tasks.Task"]) -> str:
    columns = [
        "ID", "Simulator", "Status", "Submitted", "Started", "Computation Time",
        "Resource Type"
    ]
    rows = []

    for task in tasks:
        info = task.get_info()
        simulator = task.get_simulator_name()
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
            simulator,
            status,
            format_utils.datetime_formatter(info.get("input_submit_time",
                                                     None)),
            format_utils.datetime_formatter(info.get("start_time", None)),
            execution_time,
            resource_type,
        ]
        rows.append(row)

    return columns, rows


def list_tasks(args):
    """List tasks."""

    columns, rows = _data_of_tasks(tasks.list(last_n=args.last_n))

    print(utils.format_utils.get_tasks_str(columns, rows))


def register_tasks_cli(parser):
    _cli.utils.show_help_msg(parser)

    subparsers = parser.add_subparsers()

    list_subparser = subparsers.add_parser("list", help="List tasks")
    list_subparser.add_argument(
        "-n",
        "--last-n",
        type=int,
        default=5,
        help="List last N tasks. Default: %(default)s",
    )
    # Register function to call when this subcommand is used
    list_subparser.set_defaults(func=list_tasks)
