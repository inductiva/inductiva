"""Top command for task."""
from typing import TextIO
import argparse
import sys

from inductiva import tasks


def last_modifed_file(args: argparse.Namespace, fout: TextIO = sys.stdout):
    """
    Display the last modified file for a given task.

    This function retrieves and streams information about the most recently 
    modified file associated with a specified task. It validates that the 
    task computation has started before proceeding. If the task is invalid 
    or not started, an error message is printed to `stderr`.
    """
    task_id = args.id
    task = tasks.Task(task_id)

    result, return_code = task.last_modified_file()

    if result:
        print(result, file=fout)

    return return_code


def register(parser):
    """Register the last-modifed-file command."""
    subparser = parser.add_parser("last-modified-file",
                                  help="Displays the last modified file.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva tasks last-modified-file` command displays the last"
        "modified file for the simulation. It also shows the last time the file"
        "was modified and how long it has been since that.")

    subparser.add_argument("id",
                           type=str,
                           help="The ID of the task to get the last modifed "
                           "file.")

    subparser.set_defaults(func=last_modifed_file)
