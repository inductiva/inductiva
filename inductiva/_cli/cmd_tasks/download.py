"""Download the outputs of a task by ID via CLI."""
import pathlib
import argparse

import inductiva


def download_outputs(args):
    """Download the outputs of a task by ID."""
    task_ids = args.id
    filenames = args.filenames
    output_dir = args.output_dir

    if not task_ids:
        print("No ID(s) specified.\n"
              "> Use `inductiva tasks download -h` for help.")
        return 1

    task_ids = set(task_ids)

    for task_id in task_ids:
        if output_dir is not None:
            output_dir_task = pathlib.Path(output_dir) / task_id
        else:
            output_dir_task = None
        try:
            inductiva.tasks.Task(task_id).download_outputs(
                filenames=filenames, output_dir=output_dir_task)
        except inductiva.client.exceptions.ApiException as exc:
            print(f"Error for task {task_id}:", exc)
    return 0


def register(parser):
    """Register the download task command."""

    subparser = parser.add_parser("download",
                                  help="Download tasks.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva tasks download` command allows to download the outputs"
        " of specified tasks "
        "on the platform.\n"
        "You can download multiple tasks by passing multiple ids.\n")

    subparser.add_argument("id",
                           type=str,
                           help="ID(s) of the task(s) to download.",
                           nargs="*")
    subparser.add_argument("--filenames",
                           type=str,
                           help="Names of the files to download.",
                           nargs="*")
    subparser.add_argument("--output_dir",
                           type=str,
                           help="Path of where to download the task outputs.")

    subparser.set_defaults(func=download_outputs)
