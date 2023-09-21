"""Functions to manage or retrieve user resources."""
import inductiva
import inductiva.client
from inductiva.client.apis.tags import instance_api
from inductiva.utils import format_utils


def _machine_group_list_to_str(machine_group_list) -> str:
    """Returns a string representation of a list of machine groups."""
    columns = [
        "Name", "VM Type", "# machines", "Disk Size in GB", "Spot", "Started at"
    ]
    rows = []

    for machine_group in machine_group_list:
        rows.append([
            machine_group.name, machine_group.machine_type,
            machine_group.num_machines, machine_group.disk_size_gb,
            machine_group.spot, machine_group.create_time
        ])

    formatters = {"Started at": format_utils.datetime_formatter}
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


def _fetch_machine_groups_from_api():
    """Get all active machine groups of a user from the API."""
    try:
        api = instance_api.InstanceApi(inductiva.api.get_client())
        response = api.list_active_user_instance_groups()
        if len(response.body) == 0:
            print("No active machine groups found.")
            return response.body

        return response.body

    except inductiva.client.ApiException as api_exception:
        raise api_exception


# pylint: disable=redefined-builtin
def list():
    # pylint: disable=line-too-long
    """Lists all active machine groups info.

    This method queries Google Cloud for active machine groups
    that belong to the user. It outputs all the information relative to
    each machine group as folllows:

    INFO:absl:Active machine groups:
                                        Name         VM Type   # machines    Disk Size in GB       Spot         Started at
    api-1b1f724c-5cfe-4d87-8439-9689aa139723   c2-standard-4            1                 40      False   13 Sep, 07:38:50
    api-8e6bf7d8-4888-4de9-bda5-268484b46e6f   c2-standard-4            1                 40      False   13 Sep, 07:37:49

    The name of the machine group can be used to retrieve a MachineGroup object
    with the 'get' function."""
    # pylint: enable=line-too-long

    machine_group_list = get()
    if len(machine_group_list) != 0:
        print("Active machine groups:")
        print(_machine_group_list_to_str(machine_group_list))


def get():
    """Returns a list of 'MachineGroup' objects."""

    # Retrive the active machine group names
    machine_groups = _fetch_machine_groups_from_api()

    machine_group_list = [
        inductiva.resources.MachineGroup.from_api_response(mg) \
            for mg in machine_groups
    ]

    return machine_group_list
