"""List the available simulators via CLI."""
import argparse

from inductiva.simulators import methods
from inductiva.utils import format_utils


def list_versions(args):
    """List available simulators and associated versions for each branch."""

    available = methods.list_available_images()

    if not args.dev:
        available.pop("development")
    fmtr = format_utils.get_ansi_formatter()

    def header_fmtr(s):
        return fmtr(s.upper(), format_utils.Emphasis.BOLD)

    for branch, simulators in available.items():
        sims = sorted(simulators.keys())

        data = {
            "simulator": sims,
            "versions": [", ".join(simulators[s]) for s in sims],
        }

        table = format_utils.get_tabular_str(data,
                                             header_formatters=[header_fmtr])
        header = f"Available simulators and versions for {branch} runs:"

        print(header_fmtr(header))
        print(table)


def register(parser):
    """Register the "simulators list" sub-command."""

    subparser = parser.add_parser("list",
                                  aliases=["ls"],
                                  help="List available simulators.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `list` sub-command lists all available simulators and\n"
        "associated versions for use with the Inductiva API,\n"
        "including those available for development purposes.\n")
    subparser.add_argument(
        "--dev",
        action="store_true",
        help="include development versions of the simulators.")

    subparser.set_defaults(func=list_versions)
