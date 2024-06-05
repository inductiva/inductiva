import typing_extensions

from inductiva.client.paths import PathValues
from inductiva.client.apis.paths.tasks_auth import TasksAuth
from inductiva.client.apis.paths.tasks_submit import TasksSubmit
from inductiva.client.apis.paths.tasks_task_id_input_upload_url import TasksTaskIdInputUploadUrl
from inductiva.client.apis.paths.tasks_task_id_input_uploaded import TasksTaskIdInputUploaded
from inductiva.client.apis.paths.tasks_task_id_input import TasksTaskIdInput
from inductiva.client.apis.paths.tasks_task_id import TasksTaskId
from inductiva.client.apis.paths.tasks import Tasks
from inductiva.client.apis.paths.tasks_task_id_status import TasksTaskIdStatus
from inductiva.client.apis.paths.tasks_task_id_output_list import TasksTaskIdOutputList
from inductiva.client.apis.paths.tasks_task_id_download_output_url import TasksTaskIdDownloadOutputUrl
from inductiva.client.apis.paths.tasks_task_id_output import TasksTaskIdOutput
from inductiva.client.apis.paths.tasks_task_id_resubmit import TasksTaskIdResubmit
from inductiva.client.apis.paths.tasks_task_id_kill import TasksTaskIdKill
from inductiva.client.apis.paths.tasks_task_id_disable_logs import TasksTaskIdDisableLogs
from inductiva.client.apis.paths.admin_users import AdminUsers
from inductiva.client.apis.paths.admin_users_email_api_key import AdminUsersEmailApiKey
from inductiva.client.apis.paths.admin_users_email import AdminUsersEmail
from inductiva.client.apis.paths.admin_users_email_expiry_ts import AdminUsersEmailExpiryTs
from inductiva.client.apis.paths.admin_users_username_tasks import AdminUsersUsernameTasks
from inductiva.client.apis.paths.admin_groups import AdminGroups
from inductiva.client.apis.paths.admin_groups_active import AdminGroupsActive
from inductiva.client.apis.paths.admin_groups_default import AdminGroupsDefault
from inductiva.client.apis.paths.admin_groups_default_machine_group_id import AdminGroupsDefaultMachineGroupId
from inductiva.client.apis.paths.admin_providers import AdminProviders
from inductiva.client.apis.paths.admin_providers_provider_id import AdminProvidersProviderId
from inductiva.client.apis.paths.admin_active_tasks import AdminActiveTasks
from inductiva.client.apis.paths.admin_users_username_cost import AdminUsersUsernameCost
from inductiva.client.apis.paths.admin_executer_tracker_token import AdminExecuterTrackerToken
from inductiva.client.apis.paths.executer_tracker_register import ExecuterTrackerRegister
from inductiva.client.apis.paths.executer_tracker_machine_id import ExecuterTrackerMachineId
from inductiva.client.apis.paths.executer_tracker_machine_id_task import ExecuterTrackerMachineIdTask
from inductiva.client.apis.paths.executer_tracker_machine_id_task_task_id_ack import ExecuterTrackerMachineIdTaskTaskIdAck
from inductiva.client.apis.paths.executer_tracker_machine_id_task_task_id_message import ExecuterTrackerMachineIdTaskTaskIdMessage
from inductiva.client.apis.paths.executer_tracker_machine_id_task_task_id_message_unblock import ExecuterTrackerMachineIdTaskTaskIdMessageUnblock
from inductiva.client.apis.paths.executer_tracker_machine_id_event import ExecuterTrackerMachineIdEvent
from inductiva.client.apis.paths.executer_tracker_machine_id_task_task_id_download_input_url import ExecuterTrackerMachineIdTaskTaskIdDownloadInputUrl
from inductiva.client.apis.paths.executer_tracker_machine_id_task_task_id_upload_output_url import ExecuterTrackerMachineIdTaskTaskIdUploadOutputUrl
from inductiva.client.apis.paths.executer_tracker_machine_id_task_task_id_metric import ExecuterTrackerMachineIdTaskTaskIdMetric
from inductiva.client.apis.paths.compute_group import ComputeGroup
from inductiva.client.apis.paths.compute_type import ComputeType
from inductiva.client.apis.paths.compute_group_start import ComputeGroupStart
from inductiva.client.apis.paths.compute_group_mg_id import ComputeGroupMgId
from inductiva.client.apis.paths.compute_price import ComputePrice
from inductiva.client.apis.paths.compute_groups import ComputeGroups
from inductiva.client.apis.paths.compute_group_status import ComputeGroupStatus
from inductiva.client.apis.paths.compute_machine_types import ComputeMachineTypes
from inductiva.client.apis.paths.compute_group_name import ComputeGroupName
from inductiva.client.apis.paths.storage_size import StorageSize
from inductiva.client.apis.paths.storage_contents import StorageContents
from inductiva.client.apis.paths.version import Version
from inductiva.client.apis.paths.version_check import VersionCheck
from inductiva.client.apis.paths.users_cost import UsersCost
from inductiva.client.apis.paths.users_quotas import UsersQuotas
from inductiva.client.apis.paths.users_me import UsersMe
from inductiva.client.apis.paths.projects import Projects
from inductiva.client.apis.paths.projects_name import ProjectsName

PathToApi = typing_extensions.TypedDict(
    'PathToApi', {
        PathValues.TASKS_AUTH:
            TasksAuth,
        PathValues.TASKS_SUBMIT:
            TasksSubmit,
        PathValues.TASKS_TASK_ID_INPUT_UPLOAD_URL:
            TasksTaskIdInputUploadUrl,
        PathValues.TASKS_TASK_ID_INPUT_UPLOADED:
            TasksTaskIdInputUploaded,
        PathValues.TASKS_TASK_ID_INPUT:
            TasksTaskIdInput,
        PathValues.TASKS_TASK_ID:
            TasksTaskId,
        PathValues.TASKS:
            Tasks,
        PathValues.TASKS_TASK_ID_STATUS:
            TasksTaskIdStatus,
        PathValues.TASKS_TASK_ID_OUTPUT_LIST:
            TasksTaskIdOutputList,
        PathValues.TASKS_TASK_ID_DOWNLOAD_OUTPUT_URL:
            TasksTaskIdDownloadOutputUrl,
        PathValues.TASKS_TASK_ID_OUTPUT:
            TasksTaskIdOutput,
        PathValues.TASKS_TASK_ID_RESUBMIT:
            TasksTaskIdResubmit,
        PathValues.TASKS_TASK_ID_KILL:
            TasksTaskIdKill,
        PathValues.TASKS_TASK_ID_DISABLE_LOGS:
            TasksTaskIdDisableLogs,
        PathValues.ADMIN_USERS:
            AdminUsers,
        PathValues.ADMIN_USERS_EMAIL_API_KEY:
            AdminUsersEmailApiKey,
        PathValues.ADMIN_USERS_EMAIL:
            AdminUsersEmail,
        PathValues.ADMIN_USERS_EMAIL_EXPIRY_TS:
            AdminUsersEmailExpiryTs,
        PathValues.ADMIN_USERS_USERNAME_TASKS:
            AdminUsersUsernameTasks,
        PathValues.ADMIN_GROUPS:
            AdminGroups,
        PathValues.ADMIN_GROUPS_ACTIVE:
            AdminGroupsActive,
        PathValues.ADMIN_GROUPS_DEFAULT:
            AdminGroupsDefault,
        PathValues.ADMIN_GROUPS_DEFAULT_MACHINE_GROUP_ID:
            AdminGroupsDefaultMachineGroupId,
        PathValues.ADMIN_PROVIDERS:
            AdminProviders,
        PathValues.ADMIN_PROVIDERS_PROVIDER_ID:
            AdminProvidersProviderId,
        PathValues.ADMIN_ACTIVE_TASKS:
            AdminActiveTasks,
        PathValues.ADMIN_USERS_USERNAME_COST:
            AdminUsersUsernameCost,
        PathValues.ADMIN_EXECUTERTRACKER_TOKEN:
            AdminExecuterTrackerToken,
        PathValues.EXECUTERTRACKER_REGISTER:
            ExecuterTrackerRegister,
        PathValues.EXECUTERTRACKER_MACHINE_ID:
            ExecuterTrackerMachineId,
        PathValues.EXECUTERTRACKER_MACHINE_ID_TASK:
            ExecuterTrackerMachineIdTask,
        PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_ACK:
            ExecuterTrackerMachineIdTaskTaskIdAck,
        PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_MESSAGE:
            ExecuterTrackerMachineIdTaskTaskIdMessage,
        PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_MESSAGE_UNBLOCK:
            ExecuterTrackerMachineIdTaskTaskIdMessageUnblock,
        PathValues.EXECUTERTRACKER_MACHINE_ID_EVENT:
            ExecuterTrackerMachineIdEvent,
        PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_DOWNLOAD_INPUT_URL:
            ExecuterTrackerMachineIdTaskTaskIdDownloadInputUrl,
        PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_UPLOAD_OUTPUT_URL:
            ExecuterTrackerMachineIdTaskTaskIdUploadOutputUrl,
        PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_METRIC:
            ExecuterTrackerMachineIdTaskTaskIdMetric,
        PathValues.COMPUTE_GROUP:
            ComputeGroup,
        PathValues.COMPUTE_TYPE:
            ComputeType,
        PathValues.COMPUTE_GROUP_START:
            ComputeGroupStart,
        PathValues.COMPUTE_GROUP_MG_ID:
            ComputeGroupMgId,
        PathValues.COMPUTE_PRICE:
            ComputePrice,
        PathValues.COMPUTE_GROUPS:
            ComputeGroups,
        PathValues.COMPUTE_GROUP_STATUS:
            ComputeGroupStatus,
        PathValues.COMPUTE_MACHINE_TYPES:
            ComputeMachineTypes,
        PathValues.COMPUTE_GROUP_NAME:
            ComputeGroupName,
        PathValues.STORAGE_SIZE:
            StorageSize,
        PathValues.STORAGE_CONTENTS:
            StorageContents,
        PathValues.VERSION:
            Version,
        PathValues.VERSIONCHECK:
            VersionCheck,
        PathValues.USERS_COST:
            UsersCost,
        PathValues.USERS_QUOTAS:
            UsersQuotas,
        PathValues.USERS_ME:
            UsersMe,
        PathValues.PROJECTS:
            Projects,
        PathValues.PROJECTS_NAME:
            ProjectsName,
    })

path_to_api = PathToApi({
    PathValues.TASKS_AUTH:
        TasksAuth,
    PathValues.TASKS_SUBMIT:
        TasksSubmit,
    PathValues.TASKS_TASK_ID_INPUT_UPLOAD_URL:
        TasksTaskIdInputUploadUrl,
    PathValues.TASKS_TASK_ID_INPUT_UPLOADED:
        TasksTaskIdInputUploaded,
    PathValues.TASKS_TASK_ID_INPUT:
        TasksTaskIdInput,
    PathValues.TASKS_TASK_ID:
        TasksTaskId,
    PathValues.TASKS:
        Tasks,
    PathValues.TASKS_TASK_ID_STATUS:
        TasksTaskIdStatus,
    PathValues.TASKS_TASK_ID_OUTPUT_LIST:
        TasksTaskIdOutputList,
    PathValues.TASKS_TASK_ID_DOWNLOAD_OUTPUT_URL:
        TasksTaskIdDownloadOutputUrl,
    PathValues.TASKS_TASK_ID_OUTPUT:
        TasksTaskIdOutput,
    PathValues.TASKS_TASK_ID_RESUBMIT:
        TasksTaskIdResubmit,
    PathValues.TASKS_TASK_ID_KILL:
        TasksTaskIdKill,
    PathValues.TASKS_TASK_ID_DISABLE_LOGS:
        TasksTaskIdDisableLogs,
    PathValues.ADMIN_USERS:
        AdminUsers,
    PathValues.ADMIN_USERS_EMAIL_API_KEY:
        AdminUsersEmailApiKey,
    PathValues.ADMIN_USERS_EMAIL:
        AdminUsersEmail,
    PathValues.ADMIN_USERS_EMAIL_EXPIRY_TS:
        AdminUsersEmailExpiryTs,
    PathValues.ADMIN_USERS_USERNAME_TASKS:
        AdminUsersUsernameTasks,
    PathValues.ADMIN_GROUPS:
        AdminGroups,
    PathValues.ADMIN_GROUPS_ACTIVE:
        AdminGroupsActive,
    PathValues.ADMIN_GROUPS_DEFAULT:
        AdminGroupsDefault,
    PathValues.ADMIN_GROUPS_DEFAULT_MACHINE_GROUP_ID:
        AdminGroupsDefaultMachineGroupId,
    PathValues.ADMIN_PROVIDERS:
        AdminProviders,
    PathValues.ADMIN_PROVIDERS_PROVIDER_ID:
        AdminProvidersProviderId,
    PathValues.ADMIN_ACTIVE_TASKS:
        AdminActiveTasks,
    PathValues.ADMIN_USERS_USERNAME_COST:
        AdminUsersUsernameCost,
    PathValues.ADMIN_EXECUTERTRACKER_TOKEN:
        AdminExecuterTrackerToken,
    PathValues.EXECUTERTRACKER_REGISTER:
        ExecuterTrackerRegister,
    PathValues.EXECUTERTRACKER_MACHINE_ID:
        ExecuterTrackerMachineId,
    PathValues.EXECUTERTRACKER_MACHINE_ID_TASK:
        ExecuterTrackerMachineIdTask,
    PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_ACK:
        ExecuterTrackerMachineIdTaskTaskIdAck,
    PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_MESSAGE:
        ExecuterTrackerMachineIdTaskTaskIdMessage,
    PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_MESSAGE_UNBLOCK:
        ExecuterTrackerMachineIdTaskTaskIdMessageUnblock,
    PathValues.EXECUTERTRACKER_MACHINE_ID_EVENT:
        ExecuterTrackerMachineIdEvent,
    PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_DOWNLOAD_INPUT_URL:
        ExecuterTrackerMachineIdTaskTaskIdDownloadInputUrl,
    PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_UPLOAD_OUTPUT_URL:
        ExecuterTrackerMachineIdTaskTaskIdUploadOutputUrl,
    PathValues.EXECUTERTRACKER_MACHINE_ID_TASK_TASK_ID_METRIC:
        ExecuterTrackerMachineIdTaskTaskIdMetric,
    PathValues.COMPUTE_GROUP:
        ComputeGroup,
    PathValues.COMPUTE_TYPE:
        ComputeType,
    PathValues.COMPUTE_GROUP_START:
        ComputeGroupStart,
    PathValues.COMPUTE_GROUP_MG_ID:
        ComputeGroupMgId,
    PathValues.COMPUTE_PRICE:
        ComputePrice,
    PathValues.COMPUTE_GROUPS:
        ComputeGroups,
    PathValues.COMPUTE_GROUP_STATUS:
        ComputeGroupStatus,
    PathValues.COMPUTE_MACHINE_TYPES:
        ComputeMachineTypes,
    PathValues.COMPUTE_GROUP_NAME:
        ComputeGroupName,
    PathValues.STORAGE_SIZE:
        StorageSize,
    PathValues.STORAGE_CONTENTS:
        StorageContents,
    PathValues.VERSION:
        Version,
    PathValues.VERSIONCHECK:
        VersionCheck,
    PathValues.USERS_COST:
        UsersCost,
    PathValues.USERS_QUOTAS:
        UsersQuotas,
    PathValues.USERS_ME:
        UsersMe,
    PathValues.PROJECTS:
        Projects,
    PathValues.PROJECTS_NAME:
        ProjectsName,
})
