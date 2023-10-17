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
from inductiva.client.apis.paths.admin_users import AdminUsers
from inductiva.client.apis.paths.admin_users_username_tasks import AdminUsersUsernameTasks
from inductiva.client.apis.paths.admin_tasks import AdminTasks
from inductiva.client.apis.paths.admin_groups import AdminGroups
from inductiva.client.apis.paths.executers_register import ExecutersRegister
from inductiva.client.apis.paths.gcp_instances_group import GcpInstancesGroup
from inductiva.client.apis.paths.gcp_instances_group_start import GcpInstancesGroupStart
from inductiva.client.apis.paths.gcp_instances_group_elastic import GcpInstancesGroupElastic
from inductiva.client.apis.paths.gcp_instances_price import GcpInstancesPrice
from inductiva.client.apis.paths.gcp_instances_status import GcpInstancesStatus
from inductiva.client.apis.paths.gcp_instances_group_status import GcpInstancesGroupStatus
from inductiva.client.apis.paths.gcp_instances_storage_size import GcpInstancesStorageSize
from inductiva.client.apis.paths.gcp_instances_storage_contents import GcpInstancesStorageContents
from inductiva.client.apis.paths.gcp_instances_storage_dir_name import GcpInstancesStorageDirName
from inductiva.client.apis.paths.gcp_instances_groups import GcpInstancesGroups
from inductiva.client.apis.paths.gcp_instances_group_name import GcpInstancesGroupName
from inductiva.client.apis.paths.users import Users
from inductiva.client.apis.paths.tail import Tail

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
        PathValues.ADMIN_USERS: AdminUsers,
        PathValues.ADMIN_USERS_USERNAME_TASKS: AdminUsersUsernameTasks,
        PathValues.ADMIN_TASKS: AdminTasks,
        PathValues.ADMIN_GROUPS: AdminGroups,
        PathValues.EXECUTERS_REGISTER: ExecutersRegister,
        PathValues.GCP_INSTANCES_GROUP: GcpInstancesGroup,
        PathValues.GCP_INSTANCES_GROUP_START: GcpInstancesGroupStart,
        PathValues.GCP_INSTANCES_GROUP_ELASTIC: GcpInstancesGroupElastic,
        PathValues.GCP_INSTANCES_PRICE: GcpInstancesPrice,
        PathValues.GCP_INSTANCES_STATUS: GcpInstancesStatus,
        PathValues.GCP_INSTANCES_GROUP_STATUS: GcpInstancesGroupStatus,
        PathValues.GCP_INSTANCES_STORAGE_SIZE: GcpInstancesStorageSize,
        PathValues.GCP_INSTANCES_STORAGE_CONTENTS: GcpInstancesStorageContents,
        PathValues.GCP_INSTANCES_STORAGE_DIR_NAME: GcpInstancesStorageDirName,
        PathValues.GCP_INSTANCES_GROUPS: GcpInstancesGroups,
        PathValues.GCP_INSTANCES_GROUP_NAME: GcpInstancesGroupName,
        PathValues.USERS: Users,
        PathValues.TAIL: Tail,
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
    PathValues.ADMIN_USERS: AdminUsers,
    PathValues.ADMIN_USERS_USERNAME_TASKS: AdminUsersUsernameTasks,
    PathValues.ADMIN_TASKS: AdminTasks,
    PathValues.ADMIN_GROUPS: AdminGroups,
    PathValues.EXECUTERS_REGISTER: ExecutersRegister,
    PathValues.GCP_INSTANCES_GROUP: GcpInstancesGroup,
    PathValues.GCP_INSTANCES_GROUP_START: GcpInstancesGroupStart,
    PathValues.GCP_INSTANCES_GROUP_ELASTIC: GcpInstancesGroupElastic,
    PathValues.GCP_INSTANCES_PRICE: GcpInstancesPrice,
    PathValues.GCP_INSTANCES_STATUS: GcpInstancesStatus,
    PathValues.GCP_INSTANCES_GROUP_STATUS: GcpInstancesGroupStatus,
    PathValues.GCP_INSTANCES_STORAGE_SIZE: GcpInstancesStorageSize,
    PathValues.GCP_INSTANCES_STORAGE_CONTENTS: GcpInstancesStorageContents,
    PathValues.GCP_INSTANCES_STORAGE_DIR_NAME: GcpInstancesStorageDirName,
    PathValues.GCP_INSTANCES_GROUPS: GcpInstancesGroups,
    PathValues.GCP_INSTANCES_GROUP_NAME: GcpInstancesGroupName,
    PathValues.USERS: Users,
    PathValues.TAIL: Tail,
})
