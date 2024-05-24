"""Main CLI entrypoint."""
from typing import TextIO
import argparse
import sys
import os
import re
import io

import inductiva
from inductiva import _cli
from inductiva import constants, utils
from . import loader
try:
    from . import ansi_pager
except ImportError:
    ansi_pager = None


def get_main_parser():
    parser = argparse.ArgumentParser(
        prog="inductiva",
        description="CLI tool for Inductiva API.",
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
              "will be read from the INDUCTIVA_API_KEY environment variable."),
    )

    # If no subcommand is provided, print help
    _cli.utils.show_help_msg(parser)

    # Create subcommands
    subparsers = parser.add_subparsers(title="available subcomands")

    # Load all modules starting with "cmd_" as subcommands.
    loader.load_commands(subparsers,
                         os.path.dirname(__file__),
                         "inductiva._cli",
                         prefix=constants.LOADER_COMMAND_PREFIX)
    return parser


def watch(func, every, args, cmd):
    """Run the function at regular intervals and display results in a pager."""
    cmd = re.sub(r"(-w|--watch)\s*([+-]?([0-9]*[.])?[0-9]+)?\s*", "", cmd, 1)
    header = f"> every {every}s: inductiva {cmd}"

    def action(fout: TextIO = sys.stdout):
        buffer = io.StringIO()
        func(args, fout=buffer)
        fout.clear()
        fout.write(buffer.getvalue())

    with ansi_pager.PagedOutput(header) as pager:
        scheduler = utils.scheduler.StoppableScheduler(every,
                                                       action,
                                                       args=(pager,))
        scheduler.daemon = True
        scheduler.start()
        pager.run()
        scheduler.stop()
        scheduler.join(0.5)


def main():
    parser = get_main_parser()
    args = parser.parse_args()

    if _cli.utils.check_running_for_first_time():
        answer = _cli.utils.user_autocompletion_install_prompt()
        if answer:
            _cli.utils.setup_zsh_autocompletion()

    if args.api_key:
        inductiva.set_api_key(args.api_key)

    exit_code = 0
    # Call the function associated with the subcommand
    try:
        if getattr(args, "watchable", False) and args.watch is not None:
            if ansi_pager is None:
                raise ImportError("watch module is not available.")
            else:
                cmd = " ".join(sys.argv[1:])
                watch(args.func, args.watch, args, cmd)
        else:
            exit_code = args.func(args)
    except Exception as e:  # pylint: disable=broad-except
        exit_code = 1
        print(e)
    return exit_code
