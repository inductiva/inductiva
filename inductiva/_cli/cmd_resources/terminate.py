"""CLI commands to terminate computational resources."""
import sys

from inductiva import resources
from inductiva.utils import input_functions
from ...localization import translator as __


def terminate_machine_group(args):
    """Terminate one or all computational resoruces."""
    names = args.name
    all_names = args.all
    confirm = args.confirm

    if names and all_names:
        print(
            "inductiva resources terminate: error: "
            "argument name not allowed with argument --all",
            file=sys.stderr)
        sys.exit(1)

    active_machines = resources.machine_groups.get()

    if not active_machines:
        return -1

    active_machine_names = [machine.name for machine in active_machines]

    if not all_names and not (all(item in active_machine_names
                                  for item in names)):
        print("One or more resource(s) name(s) does not exist.",
              file=sys.stderr)
        sys.exit(1)

    if not confirm:
        confirm = input_functions.user_confirmation_prompt(
            all_names, names, __("user-prompt-terminate-all"),
            __("user-prompt-terminate-big", len(names)),
            __("user-prompt-terminate-small"))
    if confirm:
        for machine in active_machines:
            if all_names:
                machine.terminate()
            elif machine.name in names:
                machine.terminate()


def register(parser):
    """Register the terminate command for the resources."""

    subparser = parser.add_parser("terminate", help="Terminate a resource.")

    subparser.add_argument("name",
                           type=str,
                           help="Name(s) of the resource(s) to terminate.",
                           nargs="*")
    subparser.add_argument("-y",
                           "--yes",
                           action="store_true",
                           dest="confirm",
                           default=False,
                           help="Sets any confirmation values to \"yes\" "
                           "automatically. Users will not be asked for "
                           "confirmation to terminate resource(s).")
    subparser.add_argument("--all",
                           action="store_true",
                           help="Terminate all machines.")

    subparser.set_defaults(func=terminate_machine_group)
