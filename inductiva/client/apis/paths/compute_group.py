from inductiva.client.paths.compute_group.post import ApiForpost
from inductiva.client.paths.compute_group.delete import ApiFordelete


class ComputeGroup(
        ApiForpost,
        ApiFordelete,
):
    pass
