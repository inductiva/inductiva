"""Download the outputs of a task by ID via CLI."""
import argparse

import inductiva


def _download(download_func, filenames):
    download_func(filenames=filenames)


def download(args):
    """Download the outputs of a task by ID."""
    task_ids = args.task_id
    filenames = args.filenames
    download_dir = args.dir

    if not args.output and not args.input:
        args.output = True

    if download_dir is not None:
        inductiva.set_output_dir(download_dir)

    if download_dir == "":
        print("`dir` must not be an empty string.")
        return 1

    if not task_ids:
        print("No ID(s) specified.\n"
              "> Use `inductiva tasks download -h` for help.")
        return 1

    task_ids = set(task_ids)

    for task_id in task_ids:
        task = inductiva.tasks.Task(task_id)
        if args.output:
            _download(task.download_outputs, filenames)
        if args.input:
            _download(task.download_inputs, filenames)
    return 0


def register(parser):
    """Register the download task command."""

    subparser = parser.add_parser("download",
                                  help="Download tasks.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "Download the input/output files of tasks with the given ID(s).\n"
        "All files are downloaded unless the --filenames list is given,\n"
        "in which case only the indicated files for each task are downloaded.\n"
        "The name of the input/output folder can be configured through the\n"
        "--dir option. Files are downloaded to a subdirectory named after the\n"
        "ID of the corresponding task. The output files are downloaded when\n"
        "the --output/-o option is passed. Similarly, the input files are\n"
        "downloaded when the --input/-i option is passed. If neither option \n"
        "is specified, the output files are downloaded by default.")

    subparser.add_argument("task_id",
                           type=str,
                           help="ID(s) of the task(s) to download.",
                           nargs="*")
    subparser.add_argument("--filenames",
                           type=str,
                           help="Names of the files to download.",
                           nargs="*")
    subparser.add_argument("--dir",
                           type=str,
                           help="Path of where to download the task"
                           "input/output files.")
    subparser.add_argument("--input",
                           "-i",
                           action="store_true",
                           help="Option to download input files.")
    subparser.add_argument("--output",
                           "-o",
                           action="store_true",
                           help="Option to download output files (default).")

    subparser.set_defaults(func=download)
