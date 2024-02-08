"""Available machine types and their number of cores."""

import os
import yaml

MACHINE_TYPES_FILE = os.path.join(os.path.dirname(__file__),
                                  "machine_types.yaml")

def get_available():
    """Gets the available machine types.
    
    Currently, this fetches the available machine types from a yaml file
    into a dictionary that can now be used to parse the information where
    needed."""

    with open(MACHINE_TYPES_FILE, "r", encoding="utf-8") as file:
        machine_types = yaml.safe_load(file)

    return machine_types


def list_available_machines():
    """List all available machines types."""

    resources_available = get_available()

    machine_types = []

    # Fetch the provider information
    for _, provider_resources in resources_available.items():
        # Fetch the available machine CPU series
        for cpu_series, series_info in provider_resources["cpu-series"].items():
            # Fetch the available RAM types and vCPUs info
            for ram_type, type_info in series_info["types"].items():
                vcpus = type_info["vcpus"]
                machine_types.extend([
                    f"{cpu_series}-{ram_type}-{vcpu}" for vcpu in vcpus])
                if type_info["lssd"]:
                    for vcpu in vcpus:
                        machine_types.append(
                            f"{cpu_series}-{ram_type}-{vcpu}-lssd")

    return tuple(machine_types)
