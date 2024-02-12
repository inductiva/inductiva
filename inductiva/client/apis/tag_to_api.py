import typing_extensions

from inductiva.client.apis.tags import TagValues
from inductiva.client.apis.tags.compute_api import ComputeApi
from inductiva.client.apis.tags.storage_api import StorageApi
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client.apis.tags.version_api import VersionApi

TagToApi = typing_extensions.TypedDict(
    'TagToApi', {
        TagValues.COMPUTE: ComputeApi,
        TagValues.STORAGE: StorageApi,
        TagValues.TASKS: TasksApi,
        TagValues.VERSION: VersionApi,
    })

tag_to_api = TagToApi({
    TagValues.COMPUTE: ComputeApi,
    TagValues.STORAGE: StorageApi,
    TagValues.TASKS: TasksApi,
    TagValues.VERSION: VersionApi,
})
