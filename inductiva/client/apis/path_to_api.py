import typing_extensions

from client.paths import PathValues
from client.apis.paths.task_submit import TaskSubmit
from client.apis.paths.task_task_id_input import TaskTaskIdInput
from client.apis.paths.task_task_id_status import TaskTaskIdStatus
from client.apis.paths.task_task_id_output import TaskTaskIdOutput
from client.apis.paths.task_task_id_kill import TaskTaskIdKill
from client.apis.paths.admin_user import AdminUser
from client.apis.paths.admin_user_username_tasks import AdminUserUsernameTasks
from client.apis.paths.admin_testing_username_tasks import AdminTestingUsernameTasks

PathToApi = typing_extensions.TypedDict(
    'PathToApi',
    {
        PathValues.TASK_SUBMIT: TaskSubmit,
        PathValues.TASK_TASK_ID_INPUT: TaskTaskIdInput,
        PathValues.TASK_TASK_ID_STATUS: TaskTaskIdStatus,
        PathValues.TASK_TASK_ID_OUTPUT: TaskTaskIdOutput,
        PathValues.TASK_TASK_ID_KILL: TaskTaskIdKill,
        PathValues.ADMIN_USER: AdminUser,
        PathValues.ADMIN_USER_USERNAME_TASKS: AdminUserUsernameTasks,
        PathValues.ADMIN_TESTING_USERNAME_TASKS: AdminTestingUsernameTasks,
    }
)

path_to_api = PathToApi(
    {
        PathValues.TASK_SUBMIT: TaskSubmit,
        PathValues.TASK_TASK_ID_INPUT: TaskTaskIdInput,
        PathValues.TASK_TASK_ID_STATUS: TaskTaskIdStatus,
        PathValues.TASK_TASK_ID_OUTPUT: TaskTaskIdOutput,
        PathValues.TASK_TASK_ID_KILL: TaskTaskIdKill,
        PathValues.ADMIN_USER: AdminUser,
        PathValues.ADMIN_USER_USERNAME_TASKS: AdminUserUsernameTasks,
        PathValues.ADMIN_TESTING_USERNAME_TASKS: AdminTestingUsernameTasks,
    }
)
