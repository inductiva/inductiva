"""List the tasks information via CLI."""
from collections import defaultdict
from typing import TextIO
import argparse
import sys

from inductiva import tasks, _cli, projects
from inductiva.utils import format_utils
from inductiva.client import models
from inductiva.utils import format_utils
from typing import Any, Iterable, List, Mapping


def _tasks_lists_to_dict(list_of_tasks: Iterable[tasks.Task]) -> Mapping[str, List[Any]]:
    """
    Converts an Iterable of tasks to a dictionary with all the
    relevant information for all the tasks.
        Args:
            list_of_tasks: An Iterable of tasks.
        Returns:
            A dictionary with all the relevant information for
            all the tasks. Example: { "ID": [1, 2, 3], 
            "Simulator": ["reef3d", "reef3d", "reef3d"], ... }
    """
    column_names = [
        "ID", "Simulator", "Status", "Submitted", "Started", "Computation Time",
        "Resource Type"
    ]
    table = defaultdict(list, {key: [] for key in column_names})

    for task in list_of_tasks:
        execution_time = task.get_computation_time(cached=True)

        if execution_time is not None:
            execution_time = format_utils.seconds_formatter(execution_time)
            if task.info.computation_end_time is None:
                if task.info.status in ["started", "submitted"]:
                    execution_time = f"*{execution_time}"
                else:
                    execution_time = "n/a"

        if task.info.executer is None:
            resource_type = None
        else:
            if task.info.executer.vm_type == "n/a":
                vm_type = task.info.executer.vm_name
            else:
                vm_type = task.info.executer.vm_type
            resource_type = (f"{task.info.executer.host_type} "
                             f"{vm_type}")
            if task.info.executer.n_mpi_hosts > 1:
                resource_type += f" x{task.info.executer.n_mpi_hosts}"

        table["ID"].append(task.id)
        table["Simulator"].append(task.get_simulator_name())
        table["Status"].append(task.info.status_alias)
        table["Submitted"].append(task.info.input_submit_time)
        table["Started"].append(task.info.start_time)
        table["Computation Time"].append(execution_time)
        table["Resource Type"].append(resource_type)

    return table

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

    table_dict = _tasks_lists_to_dict(task_list)

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
                             "the '-n' or '--last-n' option.\n")

    _cli.utils.add_watch_argument(subparser)

    group = subparser.add_mutually_exclusive_group()
    group.add_argument("-n",
                       "--last-n",
                       type=int,
                       default=10,
                       help="List last N tasks. Default: 10.")
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

    subparser.set_defaults(func=list_tasks)
