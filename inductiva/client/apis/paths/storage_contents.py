from inductiva.client.paths.storage_contents.get import ApiForget
from inductiva.client.paths.storage_contents.delete import ApiFordelete


class StorageContents(
        ApiForget,
        ApiFordelete,
):
    pass
