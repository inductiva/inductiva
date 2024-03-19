"""Downloads a tasks by id via CLI."""
import argparse
import inductiva


def download_task(args):
    """Downloads a task by id."""
    ids = args.id
    filenames = args.filenames
    output_dir = args.output_dir

    if not ids:
        print("No id(s) specified.\n"
              "> Use `inductiva tasks download -h` for help.")
        return 1

    ids = set(ids)

    if len(ids) > 1 and output_dir is not None:
        print("Cannot specify `output_dir` for more than one download")
        return 1

    for task_id in ids:
        try:
            inductiva.tasks.Task(task_id).download_outputs(
                filenames=filenames, output_dir=output_dir)
        except RuntimeError as exc:
            print(f"Error for task {task_id}:", exc)
    return 0


def register(parser):
    """Register the download task command."""

    subparser = parser.add_parser("download",
                                  help="Download tasks.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva tasks download` command downloads specified tasks "
        "on the platform.\n"
        "You can download multiple tasks by passive multiple ids.\n")

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
                           help="Where to download the task to.")

    subparser.set_defaults(func=download_task)
