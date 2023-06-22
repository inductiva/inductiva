from inductiva.client.paths.gcp_instances_group.post import ApiForpost
from inductiva.client.paths.gcp_instances_group.delete import ApiFordelete


class GcpInstancesGroup(
        ApiForpost,
        ApiFordelete,
):
    pass
