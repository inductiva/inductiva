import typing_extensions

from inductiva.client.paths import PathValues
from inductiva.client.apis.paths.task_submit import TaskSubmit
from inductiva.client.apis.paths.task_task_id_input import TaskTaskIdInput
from inductiva.client.apis.paths.task_task_id import TaskTaskId
from inductiva.client.apis.paths.task import Task
from inductiva.client.apis.paths.task_task_id_status import TaskTaskIdStatus
from inductiva.client.apis.paths.task_task_id_output import TaskTaskIdOutput
from inductiva.client.apis.paths.task_task_id_kill import TaskTaskIdKill
from inductiva.client.apis.paths.admin_users import AdminUsers
from inductiva.client.apis.paths.admin_users_username import AdminUsersUsername
from inductiva.client.apis.paths.admin_users_username_tasks import AdminUsersUsernameTasks
from inductiva.client.apis.paths.admin_tasks import AdminTasks

PathToApi = typing_extensions.TypedDict(
    'PathToApi', {
        PathValues.TASK_SUBMIT: TaskSubmit,
        PathValues.TASK_TASK_ID_INPUT: TaskTaskIdInput,
        PathValues.TASK_TASK_ID: TaskTaskId,
        PathValues.TASK: Task,
        PathValues.TASK_TASK_ID_STATUS: TaskTaskIdStatus,
        PathValues.TASK_TASK_ID_OUTPUT: TaskTaskIdOutput,
        PathValues.TASK_TASK_ID_KILL: TaskTaskIdKill,
        PathValues.ADMIN_USERS: AdminUsers,
        PathValues.ADMIN_USERS_USERNAME: AdminUsersUsername,
        PathValues.ADMIN_USERS_USERNAME_TASKS: AdminUsersUsernameTasks,
        PathValues.ADMIN_TASKS: AdminTasks,
    })

path_to_api = PathToApi({
    PathValues.TASK_SUBMIT: TaskSubmit,
    PathValues.TASK_TASK_ID_INPUT: TaskTaskIdInput,
    PathValues.TASK_TASK_ID: TaskTaskId,
    PathValues.TASK: Task,
    PathValues.TASK_TASK_ID_STATUS: TaskTaskIdStatus,
    PathValues.TASK_TASK_ID_OUTPUT: TaskTaskIdOutput,
    PathValues.TASK_TASK_ID_KILL: TaskTaskIdKill,
    PathValues.ADMIN_USERS: AdminUsers,
    PathValues.ADMIN_USERS_USERNAME: AdminUsersUsername,
    PathValues.ADMIN_USERS_USERNAME_TASKS: AdminUsersUsernameTasks,
    PathValues.ADMIN_TASKS: AdminTasks,
})
