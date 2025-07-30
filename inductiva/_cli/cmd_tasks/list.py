"""List the tasks information via CLI."""
import textwrap
from typing import TextIO
import argparse
import sys

from inductiva import tasks, _cli, projects
from inductiva.utils import format_utils
from inductiva.client import models


def list_tasks(args: argparse.Namespace, fout: TextIO = sys.stdout):
    _list_tasks(**vars(args), fout=fout)


def _list_tasks(project_name, last_n, task_id, all_tasks: bool, fout: TextIO,
                **_):
    """ List the user's tasks. 

    Lists based on the flags (task_id, last_n or project_name).
    It can print the last N tasks (default value is 10)
    or a single task with a specific ID
    or all tasks of a project.
    """
    if project_name and task_id:
        print(
            "inductiva tasks list: error: "
            "argument -i/--id not allowed with argument -p/--project-name",
            file=sys.stderr)
        return 1

    last_n = -1 if all_tasks else last_n

    if task_id:
        task_list = [tasks.Task(i) for i in task_id]
    elif project_name:
        print(f"Showing tasks for project: {project_name}.", file=fout)
        task_list = projects.Project(project_name).get_tasks(last_n=last_n)
    else:
        task_list = tasks.get_tasks(last_n=last_n)

    if not task_list:
        print("No tasks found.", file=fout)
        return 1

    table_dict = tasks.to_dict(task_list)

    emph_formatter = format_utils.get_ansi_formatter()

    def color_formater(status):
        if status.upper() == models.TaskStatusCode.SUCCESS.upper():
            return emph_formatter(status, format_utils.Emphasis.GREEN)
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
    print(
        "\nTo see more details about a task, "
        "use `inductiva tasks info <task_id>`.",)
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
                             "the `-n` or `--last-n` options.\n")

    _cli.utils.add_watch_argument(subparser)

    group = subparser.add_mutually_exclusive_group()
    group.add_argument("-n",
                       "--last-n",
                       type=int,
                       default=10,
                       help="List most recent tasks.")
    group.add_argument("-a",
                       "--all",
                       dest="all_tasks",
                       action="store_true",
                       help="List all tasks.")
    group.add_argument("-i",
                       "--id",
                       type=str,
                       nargs="+",
                       default=None,
                       dest="task_id",
                       help="List a task with a specific ID.")
    #project e id da erro exclusivo
    subparser.add_argument("-p",
                           "--project-name",
                           type=str,
                           default=None,
                           help="List the tasks of a project.")

    subparser.epilog = textwrap.dedent("""\
        examples:
            # List the last 4 tasks with their details:
            $ inductiva tasks list -n 4
            ID                              SIMULATOR          STATUS         SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
            jxwt0rm8s8xspdfcegtgkkana       splishsplash       started        08 Feb, 13:25:49       08 Feb, 13:26:04       *0:00:05                 c2-standard-4
            n0zcac8rmw7xydbis3m407kb4       splishsplash       started        08 Feb, 13:25:48       08 Feb, 13:26:03       *0:00:07                 c2-standard-4
            8nmpn4h99nyfpo4da9jw2405q       splishsplash       started        08 Feb, 13:25:47       08 Feb, 13:26:02       *0:00:09                 c2-standard-4
            so6i93pi74b89rndircubp3v2       splishsplash       started        08 Feb, 13:25:47       08 Feb, 13:26:02       *0:00:10                 c2-standard-4
    """)

    subparser.set_defaults(func=list_tasks)
