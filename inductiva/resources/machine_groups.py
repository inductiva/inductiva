"""Functions to manage or retrieve user resources."""
import inductiva
import inductiva.client
from inductiva.client.apis.tags import compute_api
from inductiva.utils import format_utils
from inductiva import resources


def _machine_group_list_to_str(machine_group_list) -> str:
    """Returns a string representation of a list of machine groups."""
    columns = [
        "Name",
        "VM Type",
        "Elastic",
        "# machines",
        "Disk Size in GB",
        "Spot",
        "Started at (UTC)",
    ]
    rows = []

    for machine_group in machine_group_list:
        if machine_group.is_elastic:
            num_active_machines = machine_group.num_active_machines
            max_machines = machine_group.max_machines
            num_active_machines = f"{num_active_machines}/{max_machines}"
        else:
            # TODO: retrieve the number of max. requested machines from the API
            num_active_machines = machine_group.num_machines
        rows.append([
            machine_group.name, machine_group.machine_type,
            machine_group.is_elastic, num_active_machines,
            machine_group.disk_size_gb, machine_group.spot,
            machine_group.create_time
        ])

    formatters = {"Started at (UTC)": format_utils.datetime_formatter}
    override_col_space = {
        "VM Type": 15,
        "# machines": 12,
        "Spot": 10,
        "Elastic": 10,
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
        api = compute_api.ComputeApi(inductiva.api.get_client())
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

    Active machine groups:
                                        Name         VM Type    Elastic   # machines    Disk Size in GB       Spot   Started at (UTC)
    api-6359c03d-c4f9-479f-8b11-ba1f8f55a58c   e2-standard-4      False            3                 40      False   10 Oct, 13:40:50
    api-db2046cf-a6fc-4124-926c-1a24329da5ea   e2-standard-4       True          2/4                 40      False   10 Oct, 12:43:03

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
    machine_group_list = []

    for mg in machine_groups:
        if mg["is_elastic"]:
            mg_class = resources.ElasticMachineGroup
        else:
            mg_class = resources.MachineGroup
        machine_group_list.append(mg_class.from_api_response(mg))

    return machine_group_list
