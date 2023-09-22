"""Main CLI entrypoint."""
import argparse

import inductiva
from inductiva import _cli


def main():
    parser = argparse.ArgumentParser(
        prog="inductiva",
        description="CLI tool for Inductiva API",
    )

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {inductiva.__version__}",
    )

    parser.add_argument(
        "--api-key",
        type=str,
        help=("API key to use. If not provided, it "
              "will be read from the INDUCTIVA_API_KEY environment variable"),
    )

    # If no subcommand is provided, print help
    _cli.utils.show_help_msg(parser)

    # Create subcommands
    subparsers = parser.add_subparsers()
    tasks_subparser = subparsers.add_parser(
        "tasks",
        help="View tasks information",
    )
    resources_subparser = subparsers.add_parser(
        "resources",
        help="View resources information",
    )

    # Register subcommands (e.g. list) for each subcommand (e.g. tasks)
    _cli.register_tasks_cli(tasks_subparser)
    _cli.register_resources_cli(resources_subparser)

    args = parser.parse_args()
    if args.api_key:
        inductiva.api_key = args.api_key

    # Call the function associated with the subcommand
    try:
        args.func(args)
    except Exception as e:  # pylint: disable=broad-except
        print(e)
