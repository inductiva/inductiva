from inductiva.client.paths.gcp_instances.post import ApiForpost
from inductiva.client.paths.gcp_instances.delete import ApiFordelete


class GcpInstances(
        ApiForpost,
        ApiFordelete,
):
    pass
