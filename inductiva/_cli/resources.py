"""Resources CLI subcommand."""
import inductiva
from inductiva import _cli


def list_resources(args):
    del args  # unused
    inductiva.resources.machine_groups.list()
    inductiva.resources.storage.get_space_used()


def list_resources_available(args):
    """List all available machines"""
    del args  # unused

    machines_dict = {
        "c2-standard-": [4, 8, 16, 30, 60],
        "c3-standard-": [4, 8, 22, 44, 88, 176],
        "c2d-standard-": [2, 4, 8, 16, 32, 56, 112],
        "c2d-highcpu-": [2, 4, 8, 16, 32, 56, 112],
        "e2-standard-": [2, 4, 8, 16, 32],
        "n2-standard-": [2, 4, 8, 16, 32, 48, 64, 80, 96, 128],
        "n2d-standard-": [2, 4, 8, 16, 32, 48, 64, 80, 96, 128, 224],
        "n1-standard-": [1, 2, 4, 8, 16, 32, 64, 96]
    }

    print("Available machine types\n")
    print("machine-type: [cores-available]")
    for machine_type, cores in machines_dict.items():
        cores_str = ", ".join(str(core) for core in cores)
        print(f"{machine_type}: [{cores_str}]")

    print("\n E.g. of machine: c2-standard-8\n")


def terminate_machine(args):
    """Terminate a machine with the given name."""
    machine_name = args.machine_name

    print("Terminating machine... If exists.")
    machines_list = inductiva.resources.machine_groups.get()

    for machine in machines_list:
        if machine.name == machine_name:
            machine.terminate()
            print(f"Terminated machine {machine_name}")
            return

    print(f"Machine {machine_name} not found.")


def register_resources_cli(parser):
    _cli.utils.show_help_msg(parser)
    subparsers = parser.add_subparsers()

    list_subparser = subparsers.add_parser(
        "list", help="List currently active resources")
    available_subparser = subparsers.add_parser(
        "available", help="List available machine types")
    terminate_subparser = subparsers.add_parser("terminate",
                                                help="Terminate a machine")

    terminate_subparser.add_argument("-n",
                                     "--machine-name",
                                     type=str,
                                     help="Name of the machine to terminate")

    # Register function to call when this subcommand is used
    list_subparser.set_defaults(func=list_resources)
    available_subparser.set_defaults(func=list_resources_available)
    terminate_subparser.set_defaults(func=terminate_machine)
