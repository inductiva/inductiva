"""Main CLI entrypoint."""
import argparse
import os

import inductiva
from inductiva import _cli
from . import loader


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
    subparsers = parser.add_subparsers(title="subcommands")

    # Load all modules starting with "cmd_" as subcommands.
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         "inductiva._cli",
                         prefix="cmd_")

    args = parser.parse_args()
    if args.api_key:
        inductiva.api_key = args.api_key

    # Call the function associated with the subcommand
    try:
        args.func(args)
    except Exception as e:  # pylint: disable=broad-except
        print(e)
