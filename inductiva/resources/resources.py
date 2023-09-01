"""Functions to manage or retrieve user resources."""
from absl import logging
import inductiva
from inductiva.client.apis.tags import instance_api
from inductiva import resources


def list_machine_groups():
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

        machine_group_names = []
        for machine_group in response.body:
            machine_group_names.append(machine_group["name"])
            logging.info("Name: %s; Number of machines: %s; Created at: %s",
                         machine_group["name"], machine_group["num_instances"],
                         machine_group["create_time"])
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
