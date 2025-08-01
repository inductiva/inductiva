"""Download the project tasks files via CLI."""
import concurrent.futures as cf
import argparse
import textwrap
from typing import List

import inductiva
from inductiva.tasks.task import Task


def download_wrapper(task: Task,
                     output_dir: str = None,
                     filenames: List[str] = None):
    """Download the project tasks files.
    This method is a wrapper around the download_outputs method of the Task
    class. We do this to be able to set the outputdir for each worker.
    """
    if output_dir:
        inductiva.set_output_dir(output_dir)
    task.download_outputs(filenames)


def download_projects(args):
    """Download the project tasks files.

    This function downloads the files of the tasks of a project.
    Args:
        args.project_name (str): Name of the project to download.
        args.output_dir (str): Directory to save the downloaded files 
            default(inductiva-output).
        args.files (list): List of files to download.
        args.std (bool): Flag to download the standard output and error files.
    The downloads are done in parallel using a ThreadPoolExecutor with a maximum
    of 10 workers.
    """

    project_name = args.project_name
    output_dir = args.output_dir
    files = args.files
    std = args.std

    if std:
        files = ["stdout.txt", "stderr.txt"]

    project = inductiva.projects.Project(project_name)
    project_tasks = project.get_tasks()

    futures = []

    with cf.ThreadPoolExecutor(max_workers=10) as executor:
        for task in project_tasks:
            future = executor.submit(download_wrapper,
                                     task=task,
                                     output_dir=output_dir,
                                     filenames=files)
            futures.append(future)

    cf.wait(futures)
    print("Downloads completed.")


def register(parser):
    """Register the projects download command."""

    subparser = parser.add_parser(
        "download",
        help="Downloads the tasks files of a project.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.add_argument("project_name",
                           type=str,
                           help="Name of the project to download.")

    subparser.description = (
        "The `inductiva projects download` command allows you to "
        "download the task outputs associated with a specific project.\n"
        "You can specify the files to download with the `--files` flag \n"
        "or download the standard output and error files with the `--std` flag."
        "\n"
        "If no files are specified, all files are downloaded.\n")

    subparser.add_argument("--output-dir",
                           type=str,
                           help="Directory to save the downloaded files.")

    group = subparser.add_mutually_exclusive_group()
    group.add_argument("--files",
                       nargs="+",
                       type=str,
                       help="Downloads the specified files from every task in "
                       "the project.")
    group.add_argument("--std",
                       action="store_true",
                       help="Downloads the standard output and error files.")

    subparser.epilog = textwrap.dedent("""\
        examples:
            # Download `stderr.txt` and `system_metrics.csv` from all tasks in
            # `my-swash-project`, and save them to the current directory.
            $ inductiva projects download my-swash-project \\
                    --files stderr.txt system_metrics.csv \\
                    --output-dir .
            Downloading simulation files to 9yhsu7cqcqemgzesuzueso1hn/outputs...
            Downloading simulation files to e5ysb582dwofj3bw6233bae99/outputs...
            Downloading simulation files to umpokj4e8rbx4379hsjdwfyuu/outputs...
            Downloading simulation files to 1k3rvmf3e0cycsiiu62umg089/outputs...
            Downloading simulation files to 9qj3fhe5tjt9bx8uvada9e8uy/outputs...
            Partial download completed to 9yhsu7cqcqemgzesuzueso1hn/outputs.
            Partial download completed to 1k3rvmf3e0cycsiiu62umg089/outputs.
            Partial download completed to umpokj4e8rbx4379hsjdwfyuu/outputs.
            Partial download completed to e5ysb582dwofj3bw6233bae99/outputs.
            Partial download completed to 9qj3fhe5tjt9bx8uvada9e8uy/outputs.
            Downloads completed.
    """)

    subparser.set_defaults(func=download_projects)
