import typing_extensions

from inductiva.client.paths import PathValues
from inductiva.client.apis.paths.tasks_submit import TasksSubmit
from inductiva.client.apis.paths.tasks_task_id_input import TasksTaskIdInput
from inductiva.client.apis.paths.tasks_task_id import TasksTaskId
from inductiva.client.apis.paths.tasks import Tasks
from inductiva.client.apis.paths.tasks_task_id_status import TasksTaskIdStatus
from inductiva.client.apis.paths.tasks_task_id_output_list import TasksTaskIdOutputList
from inductiva.client.apis.paths.tasks_task_id_output import TasksTaskIdOutput
from inductiva.client.apis.paths.tasks_task_id_kill import TasksTaskIdKill
from inductiva.client.apis.paths.admin_user import AdminUser
from inductiva.client.apis.paths.admin_user_username import AdminUserUsername
from inductiva.client.apis.paths.admin_user_username_tasks import AdminUserUsernameTasks
from inductiva.client.apis.paths.admin_groups import AdminGroups
from inductiva.client.apis.paths.executer_tracker_register import ExecuterTrackerRegister
from inductiva.client.apis.paths.compute_group import ComputeGroup
from inductiva.client.apis.paths.compute_group_start import ComputeGroupStart
from inductiva.client.apis.paths.compute_group_elastic import ComputeGroupElastic
from inductiva.client.apis.paths.compute_price import ComputePrice
from inductiva.client.apis.paths.compute_status import ComputeStatus
from inductiva.client.apis.paths.compute_group_status import ComputeGroupStatus
from inductiva.client.apis.paths.compute_groups import ComputeGroups
from inductiva.client.apis.paths.compute_group_name import ComputeGroupName
from inductiva.client.apis.paths.storage_size import StorageSize
from inductiva.client.apis.paths.storage_contents import StorageContents
from inductiva.client.apis.paths.storage_dir_name import StorageDirName
from inductiva.client.apis.paths.storage_file_tail import StorageFileTail
from inductiva.client.apis.paths.storage_file_head import StorageFileHead

PathToApi = typing_extensions.TypedDict(
    'PathToApi', {
        PathValues.TASKS_SUBMIT: TasksSubmit,
        PathValues.TASKS_TASK_ID_INPUT: TasksTaskIdInput,
        PathValues.TASKS_TASK_ID: TasksTaskId,
        PathValues.TASKS: Tasks,
        PathValues.TASKS_TASK_ID_STATUS: TasksTaskIdStatus,
        PathValues.TASKS_TASK_ID_OUTPUT_LIST: TasksTaskIdOutputList,
        PathValues.TASKS_TASK_ID_OUTPUT: TasksTaskIdOutput,
        PathValues.TASKS_TASK_ID_KILL: TasksTaskIdKill,
        PathValues.ADMIN_USER: AdminUser,
        PathValues.ADMIN_USER_USERNAME: AdminUserUsername,
        PathValues.ADMIN_USER_USERNAME_TASKS: AdminUserUsernameTasks,
        PathValues.ADMIN_GROUPS: AdminGroups,
        PathValues.EXECUTERTRACKER_REGISTER: ExecuterTrackerRegister,
        PathValues.COMPUTE_GROUP: ComputeGroup,
        PathValues.COMPUTE_GROUP_START: ComputeGroupStart,
        PathValues.COMPUTE_GROUP_ELASTIC: ComputeGroupElastic,
        PathValues.COMPUTE_PRICE: ComputePrice,
        PathValues.COMPUTE_STATUS: ComputeStatus,
        PathValues.COMPUTE_GROUP_STATUS: ComputeGroupStatus,
        PathValues.COMPUTE_GROUPS: ComputeGroups,
        PathValues.COMPUTE_GROUP_NAME: ComputeGroupName,
        PathValues.STORAGE_SIZE: StorageSize,
        PathValues.STORAGE_CONTENTS: StorageContents,
        PathValues.STORAGE_DIR_NAME: StorageDirName,
        PathValues.STORAGE_FILE_TAIL: StorageFileTail,
        PathValues.STORAGE_FILE_HEAD: StorageFileHead,
    })

path_to_api = PathToApi({
    PathValues.TASKS_SUBMIT: TasksSubmit,
    PathValues.TASKS_TASK_ID_INPUT: TasksTaskIdInput,
    PathValues.TASKS_TASK_ID: TasksTaskId,
    PathValues.TASKS: Tasks,
    PathValues.TASKS_TASK_ID_STATUS: TasksTaskIdStatus,
    PathValues.TASKS_TASK_ID_OUTPUT_LIST: TasksTaskIdOutputList,
    PathValues.TASKS_TASK_ID_OUTPUT: TasksTaskIdOutput,
    PathValues.TASKS_TASK_ID_KILL: TasksTaskIdKill,
    PathValues.ADMIN_USER: AdminUser,
    PathValues.ADMIN_USER_USERNAME: AdminUserUsername,
    PathValues.ADMIN_USER_USERNAME_TASKS: AdminUserUsernameTasks,
    PathValues.ADMIN_GROUPS: AdminGroups,
    PathValues.EXECUTERTRACKER_REGISTER: ExecuterTrackerRegister,
    PathValues.COMPUTE_GROUP: ComputeGroup,
    PathValues.COMPUTE_GROUP_START: ComputeGroupStart,
    PathValues.COMPUTE_GROUP_ELASTIC: ComputeGroupElastic,
    PathValues.COMPUTE_PRICE: ComputePrice,
    PathValues.COMPUTE_STATUS: ComputeStatus,
    PathValues.COMPUTE_GROUP_STATUS: ComputeGroupStatus,
    PathValues.COMPUTE_GROUPS: ComputeGroups,
    PathValues.COMPUTE_GROUP_NAME: ComputeGroupName,
    PathValues.STORAGE_SIZE: StorageSize,
    PathValues.STORAGE_CONTENTS: StorageContents,
    PathValues.STORAGE_DIR_NAME: StorageDirName,
    PathValues.STORAGE_FILE_TAIL: StorageFileTail,
    PathValues.STORAGE_FILE_HEAD: StorageFileHead,
})
