from inductiva.client.paths.projects.get import ApiForget
from inductiva.client.paths.projects.post import ApiForpost
from inductiva.client.paths.projects.delete import ApiFordelete


class Projects(
        ApiForget,
        ApiForpost,
        ApiFordelete,
):
    pass
