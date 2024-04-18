"""Download the outputs of a task by ID via CLI."""
import argparse

import inductiva


def download_outputs(args):
    """Download the outputs of a task by ID."""
    task_ids = args.task_id
    filenames = args.filenames
    output_dir = args.output_dir

    if output_dir is not None:
        inductiva.set_output_dir(output_dir)

    if output_dir == "":
        print("`ouput_dir` must not be an empty string.")
        return 1

    if not task_ids:
        print("No ID(s) specified.\n"
              "> Use `inductiva tasks download -h` for help.")
        return 1

    task_ids = set(task_ids)

    for task_id in task_ids:
        try:
            task = inductiva.tasks.Task(task_id)
            task.download_outputs(filenames=filenames)
        except inductiva.client.exceptions.ApiException as exc:
            print(f"Error for task {task_id}:", exc)
    return 0


def register(parser):
    """Register the download task command."""

    subparser = parser.add_parser("download",
                                  help="Download tasks.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "Download the output files produced in the context of the tasks with\n"
        "the given ID(s). All files are downloaded unless the --filenames\n"
        "list is given, in which case only the indicated files for each task\n"
        " are downloaded."
        " The name of the output folder can be configured through the\n"
        "--output_dir option. Files are downloaded to a subdirectory\n"
        "named after the ID of the corresponding task.")

    subparser.add_argument("task_id",
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
