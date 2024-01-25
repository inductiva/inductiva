"""Tasks CLI subcommands."""
from inductiva import tasks, utils
from inductiva import _cli
from inductiva.utils import format_utils


def list_tasks(args):
    """
    List tasks based on the flags (task_id, last_n).
    It can print the last N tasks (default value is 5) or a single task with a specific ID
    """
    if args.task_id is not None:
        task_list = [tasks.Task(args.task_id)]

    else:
        last_n = 5 if args.last_n is None else args.last_n
        task_list = tasks.get(last_n=last_n)

    table_dict = tasks.to_dict(task_list)

    formatters = {
        "Submitted": format_utils.datetime_formatter,
        "Started": format_utils.datetime_formatter
    }

    print(utils.format_utils.get_tabular_str(
        table_dict,
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

    list_subparser.set_defaults(func=list_tasks)
