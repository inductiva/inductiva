import typing_extensions

from inductiva.client.paths import PathValues
from inductiva.client.apis.paths.tasks_auth import TasksAuth
from inductiva.client.apis.paths.tasks_submit import TasksSubmit
from inductiva.client.apis.paths.tasks_task_id_input import TasksTaskIdInput
from inductiva.client.apis.paths.tasks_task_id import TasksTaskId
from inductiva.client.apis.paths.tasks import Tasks
from inductiva.client.apis.paths.tasks_task_id_status import TasksTaskIdStatus
from inductiva.client.apis.paths.tasks_task_id_output_list import TasksTaskIdOutputList
from inductiva.client.apis.paths.tasks_task_id_output import TasksTaskIdOutput
from inductiva.client.apis.paths.tasks_task_id_kill import TasksTaskIdKill
from inductiva.client.apis.paths.admin_users import AdminUsers
from inductiva.client.apis.paths.admin_users_email_api_key import AdminUsersEmailApiKey
from inductiva.client.apis.paths.admin_users_username_tasks import AdminUsersUsernameTasks
from inductiva.client.apis.paths.admin_groups import AdminGroups
from inductiva.client.apis.paths.admin_active_tasks import AdminActiveTasks
from inductiva.client.apis.paths.executer_tracker_register import ExecuterTrackerRegister
from inductiva.client.apis.paths.compute_group import ComputeGroup
from inductiva.client.apis.paths.compute_type import ComputeType
from inductiva.client.apis.paths.compute_group_start import ComputeGroupStart
from inductiva.client.apis.paths.compute_price import ComputePrice
from inductiva.client.apis.paths.compute_groups import ComputeGroups
from inductiva.client.apis.paths.compute_group_status import ComputeGroupStatus
from inductiva.client.apis.paths.storage_size import StorageSize
from inductiva.client.apis.paths.storage_contents import StorageContents
from inductiva.client.apis.paths.version import Version
from inductiva.client.apis.paths.version_check import VersionCheck

PathToApi = typing_extensions.TypedDict(
    'PathToApi', {
        PathValues.TASKS_AUTH: TasksAuth,
        PathValues.TASKS_SUBMIT: TasksSubmit,
        PathValues.TASKS_TASK_ID_INPUT: TasksTaskIdInput,
        PathValues.TASKS_TASK_ID: TasksTaskId,
        PathValues.TASKS: Tasks,
        PathValues.TASKS_TASK_ID_STATUS: TasksTaskIdStatus,
        PathValues.TASKS_TASK_ID_OUTPUT_LIST: TasksTaskIdOutputList,
        PathValues.TASKS_TASK_ID_OUTPUT: TasksTaskIdOutput,
        PathValues.TASKS_TASK_ID_KILL: TasksTaskIdKill,
        PathValues.ADMIN_USERS: AdminUsers,
        PathValues.ADMIN_USERS_EMAIL_API_KEY: AdminUsersEmailApiKey,
        PathValues.ADMIN_USERS_USERNAME_TASKS: AdminUsersUsernameTasks,
        PathValues.ADMIN_GROUPS: AdminGroups,
        PathValues.ADMIN_ACTIVE_TASKS: AdminActiveTasks,
        PathValues.EXECUTERTRACKER_REGISTER: ExecuterTrackerRegister,
        PathValues.COMPUTE_GROUP: ComputeGroup,
        PathValues.COMPUTE_TYPE: ComputeType,
        PathValues.COMPUTE_GROUP_START: ComputeGroupStart,
        PathValues.COMPUTE_PRICE: ComputePrice,
        PathValues.COMPUTE_GROUPS: ComputeGroups,
        PathValues.COMPUTE_GROUP_STATUS: ComputeGroupStatus,
        PathValues.STORAGE_SIZE: StorageSize,
        PathValues.STORAGE_CONTENTS: StorageContents,
        PathValues.VERSION: Version,
        PathValues.VERSIONCHECK: VersionCheck,
    })

path_to_api = PathToApi({
    PathValues.TASKS_AUTH: TasksAuth,
    PathValues.TASKS_SUBMIT: TasksSubmit,
    PathValues.TASKS_TASK_ID_INPUT: TasksTaskIdInput,
    PathValues.TASKS_TASK_ID: TasksTaskId,
    PathValues.TASKS: Tasks,
    PathValues.TASKS_TASK_ID_STATUS: TasksTaskIdStatus,
    PathValues.TASKS_TASK_ID_OUTPUT_LIST: TasksTaskIdOutputList,
    PathValues.TASKS_TASK_ID_OUTPUT: TasksTaskIdOutput,
    PathValues.TASKS_TASK_ID_KILL: TasksTaskIdKill,
    PathValues.ADMIN_USERS: AdminUsers,
    PathValues.ADMIN_USERS_EMAIL_API_KEY: AdminUsersEmailApiKey,
    PathValues.ADMIN_USERS_USERNAME_TASKS: AdminUsersUsernameTasks,
    PathValues.ADMIN_GROUPS: AdminGroups,
    PathValues.ADMIN_ACTIVE_TASKS: AdminActiveTasks,
    PathValues.EXECUTERTRACKER_REGISTER: ExecuterTrackerRegister,
    PathValues.COMPUTE_GROUP: ComputeGroup,
    PathValues.COMPUTE_TYPE: ComputeType,
    PathValues.COMPUTE_GROUP_START: ComputeGroupStart,
    PathValues.COMPUTE_PRICE: ComputePrice,
    PathValues.COMPUTE_GROUPS: ComputeGroups,
    PathValues.COMPUTE_GROUP_STATUS: ComputeGroupStatus,
    PathValues.STORAGE_SIZE: StorageSize,
    PathValues.STORAGE_CONTENTS: StorageContents,
    PathValues.VERSION: Version,
    PathValues.VERSIONCHECK: VersionCheck,
})
