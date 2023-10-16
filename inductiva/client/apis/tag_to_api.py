import typing_extensions

from inductiva.client.apis.tags import TagValues
from inductiva.client.apis.tags.admin_api import AdminApi
from inductiva.client.apis.tags.default_api import DefaultApi
from inductiva.client.apis.tags.executers_api import ExecutersApi
from inductiva.client.apis.tags.instance_api import InstanceApi
from inductiva.client.apis.tags.tasks_api import TasksApi
from inductiva.client.apis.tags.users_api import UsersApi

TagToApi = typing_extensions.TypedDict(
    'TagToApi', {
        TagValues.ADMIN: AdminApi,
        TagValues.DEFAULT: DefaultApi,
        TagValues.EXECUTERS: ExecutersApi,
        TagValues.INSTANCE: InstanceApi,
        TagValues.TASKS: TasksApi,
        TagValues.USERS: UsersApi,
    })

tag_to_api = TagToApi({
    TagValues.ADMIN: AdminApi,
    TagValues.DEFAULT: DefaultApi,
    TagValues.EXECUTERS: ExecutersApi,
    TagValues.INSTANCE: InstanceApi,
    TagValues.TASKS: TasksApi,
    TagValues.USERS: UsersApi,
})
