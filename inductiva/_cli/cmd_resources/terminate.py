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
        return 1

    active_machines = resources.machine_groups.get()

    if not active_machines:
        return 0

    # dict to map from name to machine
    name_to_machine = {machine.name: machine for machine in active_machines}
    active_machine_names = name_to_machine.keys()
    target_machine_names = set(
        names)  # the user can give the same name multiple times!!
    invalid_names = target_machine_names.difference(active_machine_names)

    if invalid_names:
        for name in invalid_names:
            print(f"Resource {name} does not exist.")
        print("Aborting.")
        return 1

    confirm = confirm or input_functions.user_confirmation_prompt(
        names, __("resources-prompt-terminate-all"),
        __("resources-prompt-terminate-big", len(names)),
        __("resources-prompt-terminate-small"), all_names)

    if not confirm:
        return 0

    machines_to_kill = active_machine_names if (all_names) else (
        target_machine_names)

    for name in machines_to_kill:
        name_to_machine[name].terminate()

    return 0


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
