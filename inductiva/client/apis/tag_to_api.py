import typing_extensions

from inductiva.client.apis.tags import TagValues
from inductiva.client.apis.tags.compute_api import ComputeApi
from inductiva.client.apis.tags.events_api import EventsApi
from inductiva.client.apis.tags.projects_api import ProjectsApi
from inductiva.client.apis.tags.simulators_api import SimulatorsApi
from inductiva.client.apis.tags.storage_api import StorageApi
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client.apis.tags.users_api import UsersApi
from inductiva.client.apis.tags.version_api import VersionApi

TagToApi = typing_extensions.TypedDict(
    'TagToApi', {
        TagValues.COMPUTE: ComputeApi,
        TagValues.EVENTS: EventsApi,
        TagValues.PROJECTS: ProjectsApi,
        TagValues.SIMULATORS: SimulatorsApi,
        TagValues.STORAGE: StorageApi,
        TagValues.TASKS: TasksApi,
        TagValues.USERS: UsersApi,
        TagValues.VERSION: VersionApi,
    })

tag_to_api = TagToApi({
    TagValues.COMPUTE: ComputeApi,
    TagValues.EVENTS: EventsApi,
    TagValues.PROJECTS: ProjectsApi,
    TagValues.SIMULATORS: SimulatorsApi,
    TagValues.STORAGE: StorageApi,
    TagValues.TASKS: TasksApi,
    TagValues.USERS: UsersApi,
    TagValues.VERSION: VersionApi,
})
