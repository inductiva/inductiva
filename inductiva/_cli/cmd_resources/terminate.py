"""CLI commands to terminate computational resources."""

from inductiva import resources


def terminate_machine_group(args):
    """Terminate one or all computational resoruces."""
    active_machines = resources.machine_groups.get()
    counter = 0

    if not active_machines:
        return

    # Doesn't run in case --all is not passed.
    # Confirm the termination of all machines.
    if args.terminate_all:
        prompt = input("Confirm the termination of all active machines? (y/n)")
        confirm = prompt.lower() in ["y", "ye", "yes"]
        if confirm:
            print("Terminating all active computational resources.")
        else:
            print("Aborting the termination of resources.")
            return

    for machine in active_machines:
        if args.terminate_all:
            machine.terminate()
        elif machine.name == args.name:
            machine.terminate()
            break
        else:
            counter += 1

    if counter == len(active_machines):
        print("No active computational resources found with the "
              f"name {args.name}.")


def register(parser):
    """Register the terminate command for the resources."""

    subparser = parser.add_parser("terminate",
                                  description="Terminate resources.",
                                  help="Terminate resources.")

    group = subparser.add_mutually_exclusive_group(required=True)
    group.add_argument("-n",
                       "--name",
                       type=str,
                       help="Name of the resource to terminate.")
    group.add_argument("--all",
                       action="store_true",
                       dest="terminate_all",
                       help="Terminate all machines.")

    subparser.set_defaults(func=terminate_machine_group)
