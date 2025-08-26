"""List the available simulators via CLI."""
import argparse
import textwrap

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

    subparser.description = textwrap.dedent("""\
        The `inductiva simulators list` sub-command lists all available
        simulators and associated versions for use with the Inductiva API,
        including those available for development purposes.
    """)

    subparser.add_argument(
        "--dev",
        action="store_true",
        help="include development versions of the simulators.")

    subparser.epilog = textwrap.dedent("""\
        examples:
            $ inductiva simulators list
            AVAILABLE SIMULATORS AND VERSIONS FOR PRODUCTION RUNS:

            SIMULATOR             VERSIONS
            amr-wind              1.4.0, 1.4.0_gpu, 3.4.0, 3.4.0_gpu, 3.4.1, 3.4.1_gpu
            apptainer-converter   0.1.0
            cans                  2.3.4, 2.4.0, 2.4.0_gpu, 3.0.0, 3.0.0_gpu
            cm1                   18, 21.1
            coawst                3.8
            cp2k                  2025.1, 2025.1_gpu
            delft3d               6.04.00
            dualsphysics          5.2.1, 5.2.1_gpu, 5.4.1, 5.4.1_gpu
            ...
    """)

    subparser.set_defaults(func=list_versions)
