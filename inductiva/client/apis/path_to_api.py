import typing_extensions

from inductiva.client.paths import PathValues
from inductiva.client.apis.paths.tasks_submit import TasksSubmit
from inductiva.client.apis.paths.tasks_task_id_input import TasksTaskIdInput
from inductiva.client.apis.paths.tasks_task_id import TasksTaskId
from inductiva.client.apis.paths.tasks import Tasks
from inductiva.client.apis.paths.tasks_task_id_status import TasksTaskIdStatus
from inductiva.client.apis.paths.tasks_task_id_output import TasksTaskIdOutput
from inductiva.client.apis.paths.tasks_task_id_kill import TasksTaskIdKill
from inductiva.client.apis.paths.admin_users import AdminUsers
from inductiva.client.apis.paths.admin_users_username import AdminUsersUsername
from inductiva.client.apis.paths.admin_users_username_tasks import AdminUsersUsernameTasks
from inductiva.client.apis.paths.admin_tasks import AdminTasks
from inductiva.client.apis.paths.executers_register import ExecutersRegister
from inductiva.client.apis.paths.executers_pools import ExecutersPools
from inductiva.client.apis.paths.gcp_instances import GcpInstances
from inductiva.client.apis.paths.gcp_instances_group import GcpInstancesGroup
from inductiva.client.apis.paths.gcp_instances_price import GcpInstancesPrice

PathToApi = typing_extensions.TypedDict(
    'PathToApi', {
        PathValues.TASKS_SUBMIT: TasksSubmit,
        PathValues.TASKS_TASK_ID_INPUT: TasksTaskIdInput,
        PathValues.TASKS_TASK_ID: TasksTaskId,
        PathValues.TASKS: Tasks,
        PathValues.TASKS_TASK_ID_STATUS: TasksTaskIdStatus,
        PathValues.TASKS_TASK_ID_OUTPUT: TasksTaskIdOutput,
        PathValues.TASKS_TASK_ID_KILL: TasksTaskIdKill,
        PathValues.ADMIN_USERS: AdminUsers,
        PathValues.ADMIN_USERS_USERNAME: AdminUsersUsername,
        PathValues.ADMIN_USERS_USERNAME_TASKS: AdminUsersUsernameTasks,
        PathValues.ADMIN_TASKS: AdminTasks,
        PathValues.EXECUTERS_REGISTER: ExecutersRegister,
        PathValues.EXECUTERS_POOLS: ExecutersPools,
        PathValues.GCP_INSTANCES: GcpInstances,
        PathValues.GCP_INSTANCES_GROUP: GcpInstancesGroup,
        PathValues.GCP_INSTANCES_PRICE: GcpInstancesPrice,
    })

path_to_api = PathToApi({
    PathValues.TASKS_SUBMIT: TasksSubmit,
    PathValues.TASKS_TASK_ID_INPUT: TasksTaskIdInput,
    PathValues.TASKS_TASK_ID: TasksTaskId,
    PathValues.TASKS: Tasks,
    PathValues.TASKS_TASK_ID_STATUS: TasksTaskIdStatus,
    PathValues.TASKS_TASK_ID_OUTPUT: TasksTaskIdOutput,
    PathValues.TASKS_TASK_ID_KILL: TasksTaskIdKill,
    PathValues.ADMIN_USERS: AdminUsers,
    PathValues.ADMIN_USERS_USERNAME: AdminUsersUsername,
    PathValues.ADMIN_USERS_USERNAME_TASKS: AdminUsersUsernameTasks,
    PathValues.ADMIN_TASKS: AdminTasks,
    PathValues.EXECUTERS_REGISTER: ExecutersRegister,
    PathValues.EXECUTERS_POOLS: ExecutersPools,
    PathValues.GCP_INSTANCES: GcpInstances,
    PathValues.GCP_INSTANCES_GROUP: GcpInstancesGroup,
    PathValues.GCP_INSTANCES_PRICE: GcpInstancesPrice,
})
