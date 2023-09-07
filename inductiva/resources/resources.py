"""Functions to manage or retrieve user resources."""
from absl import logging
import inductiva
from inductiva.client.apis.tags import instance_api
from inductiva import resources
from inductiva.utils import format_utils


def _machine_group_list_to_str(machine_groups) -> str:
    """Returns a string representation of a list of machine groups."""
    columns = ["Name", "VM Type", "# machines", "Created at"]
    rows = []

    for machine_group in machine_groups:
        rows.append([
            machine_group["name"], machine_group["machine_type"],
            machine_group["num_instances"], machine_group["create_time"]
        ])

    formatters = {"Created at": format_utils.datetime_formatter}
    override_col_space = {
        "VM Type": 15,
        "# machines": 12,
    }

    return format_utils.get_tabular_str(
        rows,
        columns,
        default_col_space=18,
        override_col_space=override_col_space,
        formatters=formatters,
    )


def list_active_machine_groups():
    """Lists all active machine group names.

    Returns:
        List of machine group names.
    """
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        response = api.list_active_user_instance_groups()
        if len(response.body) == 0:
            logging.info("No active machine groups found.")
            return response.body

        machine_groups = list(response.body)
        machine_group_names = [
            machine_group["name"] for machine_group in machine_groups
        ]

        logging.info("Active machine groups:\n%s",
                     _machine_group_list_to_str(machine_groups))

        return machine_group_names

    except inductiva.client.ApiException as api_exception:
        raise api_exception


def get_machine_group(name: str):
    """Returns a 'MachineGroup' object.

    Given the name of the machine group, returns a 'MachineGroup' object
    with the attributes of that machine group.

    Args:
        name: Name of the machine group."""
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        response = api.get_instance_group({"name": name})

        mg = resources.MachineGroup(machine_type=response.body["machine_type"],
                                    num_machines=response.body["num_instances"],
                                    spot=response.body["spot"],
                                    disk_size_gb=response.body["disk_size_gb"],
                                    zone=response.body["zone"])
        mg.name = response.body["name"]
        return mg

    except inductiva.client.ApiException as api_exception:
        raise api_exception
