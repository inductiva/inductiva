"""List the tasks information via CLI."""
from typing import TextIO
import argparse
import sys

from inductiva import tasks, _cli, projects
from inductiva.utils import format_utils
from inductiva.client import models


def list_tasks(args, fout: TextIO = sys.stdout):
    """ List the last user's tasks. 

    Lists based on the flags (task_id, last_n or project_name).
    It can print the last N tasks (default value is 5)
    or a single task with a specific ID
    or all tasks of a project.
    """
    project_name = args.project_name
    task_id = args.task_id
    last_n = args.last_n

    if task_id is not None and (last_n is not None or project_name is not None):
        print(
            "inductiva tasks list: error: "
            "argument id not allowed with argument --project-name or --last-n",
            file=sys.stderr)
        return 1

    if project_name is not None:
        print(f"Showing tasks for project: {project_name}.")
        # With project the default last_n is -1 (all tasks)
        last_n = -1 if args.last_n is None else args.last_n
        task_list = projects.Project(project_name).get_tasks(last_n=last_n)
    elif task_id is not None:
        task_list = [tasks.Task(task_id)]
    else:
        # With no project the default last_n is 10
        last_n = 10 if args.last_n is None else args.last_n
        task_list = tasks.get_tasks(last_n=last_n)

    if not task_list:
        print("No tasks found.", file=fout)
        return 1

    table_dict = tasks.to_dict(task_list)

    emph_formatter = format_utils.get_ansi_formatter()

    def color_formater(status):
        if status == models.TaskStatusCode.SUCCESS:
            return emph_formatter(status, format_utils.Emphasis.GREEN)
        elif status in tasks.Task.FAILED_STATUSES:
            return emph_formatter(status, format_utils.Emphasis.RED)
        return status

    formatters = {
        "Submitted": [format_utils.datetime_formatter,],
        "Started": [format_utils.datetime_formatter,],
        "Status": [color_formater]
    }

    header_formatters = [
        lambda x: emph_formatter(x.upper(), format_utils.Emphasis.BOLD)
    ]

    print(format_utils.get_tabular_str(table_dict,
                                       formatters=formatters,
                                       header_formatters=header_formatters),
          file=fout,
          end="")

    return 0


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

    _cli.utils.add_watch_argument(subparser)

    subparser.add_argument("-n",
                           "--last-n",
                           type=int,
                           help="List last N tasks. Default: 5.")
    subparser.add_argument("-id",
                           "--task-id",
                           type=str,
                           help="List a task with a specific ID.")
    subparser.add_argument("-p",
                           "--project-name",
                           type=str,
                           help="List the tasks of a project.")

    subparser.set_defaults(func=list_tasks)
