"""Download the project tasks files via CLI."""
import concurrent.futures as cf
import argparse

import inductiva


def download_projects(args):
    """Download the project tasks files.

    This function downloads the files of the tasks of a project. The project
    name is specified with the --project-name flag. The files to download can
    be specified with the --files flag. If no files are specified, all files are
    downloaded. The standard output and error files can be downloaded with the
    --std flag. The files are saved in the default directory (inductiva-output)
    if no output directory is specified with the --output-dir flag.

    The downloads are done in parallel using a ThreadPoolExecutor with a maximum
    of 100 workers.
    """

    project_name = args.project_name
    output_dir = args.output_dir
    files = args.files
    std = args.std

    if std:
        files = ["stdout.txt", "stderr.txt"]

    project = inductiva.projects.Project(project_name)
    project_tasks = project.get_tasks()

    if output_dir is not None:
        inductiva.set_output_dir(output_dir)

    futures = []

    with cf.ThreadPoolExecutor(max_workers=100) as executor:
        for task in project_tasks:
            future = executor.submit(task.download_outputs, filenames=files)
            futures.append(future)

    cf.wait(futures)
    print("Downloads completed.")


def register(parser):
    """Register the projetcs list command."""

    subparser = parser.add_parser("download",
                                  help="List the user's projetcs.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.add_argument("project_name",
                           type=str,
                           help="Name of the project to download.")

    subparser.description = (
        "The `inductiva projects download` command provides "
        "a way to download the tasks files of your projects.\n"
        "You can specify the files to download with the --files flag \n"
        "or download the standard output and error files with the --std flag.\n"
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

    subparser.set_defaults(func=download_projects)
