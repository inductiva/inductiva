"""Download the outputs of a task by ID via CLI."""
import argparse
import textwrap

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

    subparser.description = textwrap.dedent("""\
        Download the input and/or output files of tasks with the given ID(s).

        By default, all output files of the provided task are downloaded to a
        local directory at `inductiva_output/<TASK_ID>/outputs`, relative to
        your current working directory.
        
        To download specific files, use the `--filenames` option and provide a
        list of filenames.
                                            
        You can specify multiple task IDs to download files from several tasks
        at once.                                    

        The target directory can be set with the `--dir` option. Files will be
        saved into a subdirectory named after each task ID.

        Use the `--output` (`-o`) option to download output files, and the
        `--input` (`-i`) option to download input files. If neither is
        specified, only the output files are downloaded by default.
    """)

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
                           help="Path to the directory where input/output "
                           "files will be downloaded.")
    subparser.add_argument("--input",
                           "-i",
                           action="store_true",
                           help="Option to download input files.")
    subparser.add_argument("--output",
                           "-o",
                           action="store_true",
                           help="Option to download output files.")

    subparser.set_defaults(func=download)
