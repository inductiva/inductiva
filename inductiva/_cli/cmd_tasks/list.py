"""List the tasks information via CLI."""
import argparse

from inductiva import tasks, utils
from inductiva.client import models


def list_tasks(args):
    """ List the last user's tasks. 

    Lists based on the flags (task_id, last_n).
    It can print the last N tasks (default value is 5)
    or a single task with a specific ID
    """
    if args.task_id is not None:
        task_list = [tasks.Task(args.task_id)]

    else:
        last_n = 5 if args.last_n is None else args.last_n
        task_list = tasks.get(last_n=last_n)

    table_dict = tasks.to_dict(task_list)

    emph_formatter = utils.format_utils.get_ansi_formatter()

    def color_formater(status):
        if status == models.TaskStatusCode.SUCCESS:
            return emph_formatter(status, utils.format_utils.Emphasis.GREEN)
        elif status in tasks.Task.FAILED_STATUSES:
            return emph_formatter(status, utils.format_utils.Emphasis.RED)
        return status

    formatters = {
        "Submitted": [utils.format_utils.datetime_formatter,],
        "Started": [utils.format_utils.datetime_formatter,],
        "Status": [color_formater]
    }

    header_formatters = [
        lambda x: emph_formatter(x.upper(), utils.format_utils.Emphasis.BOLD)
    ]

    print(
        utils.format_utils.get_tabular_str(table_dict,
                                           formatters=formatters,
                                           header_formatters=header_formatters))


def register(parser):
    """Register the list user's tasks command."""

    subparser = parser.add_parser("list",
                                  help="List tasks.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva tasks list` command provides "
                             "an overview of your tasks on the platform.\n"
                             "It lists the most recent tasks or a specific "
                             "task by ID.\n"
                             "You can control the number of tasks listed with "
                             "the '-n' or '--last-n' option.\n")

    group = subparser.add_mutually_exclusive_group()
    group.add_argument("-n",
                       "--last-n",
                       type=int,
                       help="List last N tasks. Default: 5.")
    group.add_argument("-id",
                       "--task-id",
                       type=str,
                       help="List a task with a specific ID.")

    subparser.set_defaults(func=list_tasks)
