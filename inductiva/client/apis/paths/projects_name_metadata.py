from inductiva.client.paths.projects_name_metadata.get import ApiForget
from inductiva.client.paths.projects_name_metadata.put import ApiForput
from inductiva.client.paths.projects_name_metadata.patch import ApiForpatch


class ProjectsNameMetadata(
        ApiForget,
        ApiForput,
        ApiForpatch,
):
    pass
