import typing_extensions

from inductiva.client.apis.tags import TagValues
from inductiva.client.apis.tags.admin_api import AdminApi
from inductiva.client.apis.tags.executers_api import ExecutersApi
from inductiva.client.apis.tags.tasks_api import TasksApi

TagToApi = typing_extensions.TypedDict(
    'TagToApi', {
        TagValues.ADMIN: AdminApi,
        TagValues.EXECUTERS: ExecutersApi,
        TagValues.TASKS: TasksApi,
    })

tag_to_api = TagToApi({
    TagValues.ADMIN: AdminApi,
    TagValues.EXECUTERS: ExecutersApi,
    TagValues.TASKS: TasksApi,
})
