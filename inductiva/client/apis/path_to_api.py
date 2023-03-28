import typing_extensions

from inductiva.client.paths import PathValues
from inductiva.client.apis.paths.task_submit import TaskSubmit
from inductiva.client.apis.paths.task_task_id_input import TaskTaskIdInput
from inductiva.client.apis.paths.task_task_id_status import TaskTaskIdStatus
from inductiva.client.apis.paths.task_task_id_output import TaskTaskIdOutput
from inductiva.client.apis.paths.task_task_id_kill import TaskTaskIdKill
from inductiva.client.apis.paths.admin_user import AdminUser
from inductiva.client.apis.paths.admin_user_username_tasks import AdminUserUsernameTasks

PathToApi = typing_extensions.TypedDict(
    'PathToApi', {
        PathValues.TASK_SUBMIT: TaskSubmit,
        PathValues.TASK_TASK_ID_INPUT: TaskTaskIdInput,
        PathValues.TASK_TASK_ID_STATUS: TaskTaskIdStatus,
        PathValues.TASK_TASK_ID_OUTPUT: TaskTaskIdOutput,
        PathValues.TASK_TASK_ID_KILL: TaskTaskIdKill,
        PathValues.ADMIN_USER: AdminUser,
        PathValues.ADMIN_USER_USERNAME_TASKS: AdminUserUsernameTasks,
    })

path_to_api = PathToApi({
    PathValues.TASK_SUBMIT: TaskSubmit,
    PathValues.TASK_TASK_ID_INPUT: TaskTaskIdInput,
    PathValues.TASK_TASK_ID_STATUS: TaskTaskIdStatus,
    PathValues.TASK_TASK_ID_OUTPUT: TaskTaskIdOutput,
    PathValues.TASK_TASK_ID_KILL: TaskTaskIdKill,
    PathValues.ADMIN_USER: AdminUser,
    PathValues.ADMIN_USER_USERNAME_TASKS: AdminUserUsernameTasks,
})
