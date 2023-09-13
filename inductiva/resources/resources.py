"""Functions to manage or retrieve user resources."""
from absl import logging
import inductiva
from inductiva.client.apis.tags import instance_api
from inductiva import resources
from inductiva.utils import format_utils


def _machine_group_list_to_str(machine_groups) -> str:
    """Returns a string representation of a list of machine groups."""
    columns = [
        "Name", "VM Type", "# machines", "Disk Size in GB", "Spot", "Created at"
    ]
    rows = []

    for machine_group in machine_groups:
        rows.append([
            machine_group["name"], machine_group["machine_type"],
            machine_group["num_instances"], machine_group["disk_size_gb"],
            bool(machine_group["spot"]), machine_group["create_time"]
        ])

    formatters = {"Created at": format_utils.datetime_formatter}
    override_col_space = {
        "VM Type": 15,
        "# machines": 12,
        "Spot": 10,
    }

    return format_utils.get_tabular_str(
        rows,
        columns,
        default_col_space=18,
        override_col_space=override_col_space,
        formatters=formatters,
    )


def list_active_machine_groups():
    # pylint: disable=line-too-long
    """Lists all active machine groups info.

    This method queries Google Cloud for active machine groups
    that belong to the user. It outputs all the information relative to
    each machine group as folllows:

    INFO:absl:Active machine groups:
                                        Name         VM Type   # machines    Disk Size in GB       Spot         Created at
    api-1b1f724c-5cfe-4d87-8439-9689aa139723   c2-standard-4            1                 40      False   13 Sep, 07:38:50
    api-8e6bf7d8-4888-4de9-bda5-268484b46e6f   c2-standard-4            1                 40      False   13 Sep, 07:37:49

    The name of the machine group can be used to retrieve a MachineGroup object
    with the 'get_machine_group' function.

    Returns:
        List of machine group names.
    """
    # pylint: enable=line-too-long
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


def get_machine_groups():
    """Returns a list of 'MachineGroup' objects."""

    # Retrive the active machine group names
    machine_group_name = list_active_machine_groups()
    machine_groups = [get_machine_group(name) for name in machine_group_name]

    return machine_groups
