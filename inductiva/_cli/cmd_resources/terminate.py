"""CLI commands to terminate computational resources."""
import sys

import argparse

from inductiva import resources
from inductiva.utils import input_functions
from ...localization import translator as __

import concurrent.futures


def terminate_machine_group(args):
    """Terminate one or all computational resources."""
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
        print("No active resources to terminate.")
        return 0

    if not all_names and not names:
        print("No resource(s) specified.\n"
              "> Use `inductiva resources terminate -h` for help.")
        return 1

    # dict to map from name to machine
    name_to_machine = {machine.name: machine for machine in active_machines}
    # the user can give the same name multiple times
    target_machine_names = set(names or name_to_machine.keys())
    invalid_names = target_machine_names.difference(name_to_machine)

    if invalid_names:
        for name in invalid_names:
            print(f"Resource {name} does not exist.")
        print("Aborting.")
        return 1

    confirm = confirm or input_functions.user_confirmation_prompt(
        target_machine_names,
        __("resources-prompt-terminate-big", len(target_machine_names)),
        __("resources-prompt-terminate-small"))

    if not confirm:
        return 0

    def terminate_machine(name):
        name_to_machine[name].terminate()

    # if there is only one machine to terminate do it in the main thread
    if len(target_machine_names) == 1:
        terminate_machine(target_machine_names.pop())
        return 0

    with concurrent.futures.ThreadPoolExecutor(
            max_workers=min(32, len(target_machine_names))) as executor:
        futures = [
            executor.submit(terminate_machine, name)
            for name in target_machine_names
        ]
        concurrent.futures.wait(futures)

    return 0


def register(parser):
    """Register the terminate command for the resources."""

    subparser = parser.add_parser("terminate",
                                  help="Terminate resources.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = ("The `inductiva resources terminate` command "
                             "provides a utility for terminating\n"
                             "active computational resources. It allows you"
                             " to specify the names of the resources\n"
                             "to terminate, or terminate all active resources."
                             " Multiple resources can be terminated\n"
                             "at once by providing their names.\n\n")

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
