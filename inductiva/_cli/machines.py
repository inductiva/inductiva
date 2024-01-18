"""CLI utilities for machines."""
import inductiva


def list_machine_groups(args):
    """List machine groups."""
    del args  # unused
    inductiva.resources.machine_groups.list()


def start_machine_group(args):
    """Start a machine group."""
    machine_type = args.machine_type
    num_machines = args.num_machines
    disk_size_gb = args.disk_size
    spot = args.spot

    machine = inductiva.resources.MachineGroup(machine_type=machine_type,
                                               num_machines=num_machines,
                                               disk_size_gb=disk_size_gb,
                                               spot=spot)

    machine.start()

    print(f"Machine group started.\nName: {machine.name}")


def list_machine_types_available(args):
    """List all available machines"""
    del args  # unused

    print("Available machine types\n")
    print("machine-type: [cores-available]")
    for (machine_type,
         cores) in inductiva.resources.machine_types.AVAILABLE_MACHINES.items():
        cores_str = ", ".join(str(core) for core in cores)
        print(f"{machine_type}: [{cores_str}]")

    print("\n E.g. of machine-type: c2-standard-8\n")


def estimate_machine_cost(args):
    machine_type = args.machine_type
    spot = args.spot

    cost = inductiva.resources.estimate_machine_cost(
        machine_type=machine_type,
        spot=spot,
    )
    print(f"Estimated cost of machine: {cost} $/h.")


def terminate_machine_group(args):
    """Terminate a machine group from a given name."""
    machine_name = args.name

    print("Terminating machine group... If exists.")
    machines_list = inductiva.resources.machine_groups.get()

    for machine in machines_list:
        if machine.name == machine_name:
            machine.terminate()
            print(f"Terminated machine group: {machine_name}")
            return

    print(f"Machine {machine_name} not found.")
