import typing_extensions

from inductiva.client.paths import PathValues
from inductiva.client.apis.paths.tasks_auth import TasksAuth
from inductiva.client.apis.paths.tasks_submit import TasksSubmit
from inductiva.client.apis.paths.tasks_task_id_input_upload_url import TasksTaskIdInputUploadUrl
from inductiva.client.apis.paths.tasks_task_id_input_uploaded import TasksTaskIdInputUploaded
from inductiva.client.apis.paths.tasks_task_id import TasksTaskId
from inductiva.client.apis.paths.tasks import Tasks
from inductiva.client.apis.paths.tasks_task_id_status import TasksTaskIdStatus
from inductiva.client.apis.paths.tasks_task_id_output_list import TasksTaskIdOutputList
from inductiva.client.apis.paths.tasks_task_id_download_input_url import TasksTaskIdDownloadInputUrl
from inductiva.client.apis.paths.tasks_task_id_download_output_url import TasksTaskIdDownloadOutputUrl
from inductiva.client.apis.paths.tasks_task_id_resubmit import TasksTaskIdResubmit
from inductiva.client.apis.paths.tasks_task_id_kill import TasksTaskIdKill
from inductiva.client.apis.paths.tasks_task_id_disable_logs import TasksTaskIdDisableLogs
from inductiva.client.apis.paths.tasks_task_id_files import TasksTaskIdFiles
from inductiva.client.apis.paths.tasks_task_id_register import TasksTaskIdRegister
from inductiva.client.apis.paths.tasks_task_id_offer import TasksTaskIdOffer
from inductiva.client.apis.paths.tasks_task_id_message import TasksTaskIdMessage
from inductiva.client.apis.paths.admin_users import AdminUsers
from inductiva.client.apis.paths.admin_users_email_terms_and_conditions import AdminUsersEmailTermsAndConditions
from inductiva.client.apis.paths.admin_users_username_organization import AdminUsersUsernameOrganization
from inductiva.client.apis.paths.admin_users_username_tier import AdminUsersUsernameTier
from inductiva.client.apis.paths.admin_users_username_credits import AdminUsersUsernameCredits
from inductiva.client.apis.paths.admin_users_email_api_key import AdminUsersEmailApiKey
from inductiva.client.apis.paths.admin_users_email import AdminUsersEmail
from inductiva.client.apis.paths.admin_users_username_storage_size_fs import AdminUsersUsernameStorageSizeFs
from inductiva.client.apis.paths.admin_users_username_storage_size import AdminUsersUsernameStorageSize
from inductiva.client.apis.paths.admin_users_username_tasks import AdminUsersUsernameTasks
from inductiva.client.apis.paths.admin_users_username_capabilities import AdminUsersUsernameCapabilities
from inductiva.client.apis.paths.admin_groups import AdminGroups
from inductiva.client.apis.paths.admin_groups_active import AdminGroupsActive
from inductiva.client.apis.paths.admin_active_tasks import AdminActiveTasks
from inductiva.client.apis.paths.admin_groups_machine_group_id_terminate import AdminGroupsMachineGroupIdTerminate
from inductiva.client.apis.paths.admin_organizations import AdminOrganizations
from inductiva.client.apis.paths.admin_organizations_organization_id import AdminOrganizationsOrganizationId
from inductiva.client.apis.paths.admin_organizations_costs import AdminOrganizationsCosts
from inductiva.client.apis.paths.admin_tiers import AdminTiers
from inductiva.client.apis.paths.admin_terminate_machine_groups_credits_exhausted import AdminTerminateMachineGroupsCreditsExhausted
from inductiva.client.apis.paths.task_runner_register import TaskRunnerRegister
from inductiva.client.apis.paths.task_runner_machine_id import TaskRunnerMachineId
from inductiva.client.apis.paths.task_runner_machine_id_task import TaskRunnerMachineIdTask
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_message import TaskRunnerMachineIdTaskTaskIdMessage
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_message_unblock import TaskRunnerMachineIdTaskTaskIdMessageUnblock
from inductiva.client.apis.paths.task_runner_machine_id_event import TaskRunnerMachineIdEvent
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_operation import TaskRunnerMachineIdTaskTaskIdOperation
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_operation_operation_id_done import TaskRunnerMachineIdTaskTaskIdOperationOperationIdDone
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_download_input_url import TaskRunnerMachineIdTaskTaskIdDownloadInputUrl
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_upload_output_url import TaskRunnerMachineIdTaskTaskIdUploadOutputUrl
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_metric import TaskRunnerMachineIdTaskTaskIdMetric
from inductiva.client.apis.paths.task_runner_machine_id_resize_disk import TaskRunnerMachineIdResizeDisk
from inductiva.client.apis.paths.task_runner_machine_id_resize_disk_done import TaskRunnerMachineIdResizeDiskDone
from inductiva.client.apis.paths.task_runner_machine_id_download_urls import TaskRunnerMachineIdDownloadUrls
from inductiva.client.apis.paths.compute_group import ComputeGroup
from inductiva.client.apis.paths.compute_group_start import ComputeGroupStart
from inductiva.client.apis.paths.compute_price import ComputePrice
from inductiva.client.apis.paths.compute_groups import ComputeGroups
from inductiva.client.apis.paths.compute_groups_history import ComputeGroupsHistory
from inductiva.client.apis.paths.compute_group_status import ComputeGroupStatus
from inductiva.client.apis.paths.compute_machine_types import ComputeMachineTypes
from inductiva.client.apis.paths.compute_group_name import ComputeGroupName
from inductiva.client.apis.paths.storage_size import StorageSize
from inductiva.client.apis.paths.storage_cost import StorageCost
from inductiva.client.apis.paths.storage_contents import StorageContents
from inductiva.client.apis.paths.storage_input_url import StorageInputUrl
from inductiva.client.apis.paths.storage_signed_urls import StorageSignedUrls
from inductiva.client.apis.paths.storage_input_notify import StorageInputNotify
from inductiva.client.apis.paths.storage_input_remote import StorageInputRemote
from inductiva.client.apis.paths.storage_ import Storage
from inductiva.client.apis.paths.storage_export import StorageExport
from inductiva.client.apis.paths.storage_operations_operation_id import StorageOperationsOperationId
from inductiva.client.apis.paths.storage_operations import StorageOperations
from inductiva.client.apis.paths.version import Version
from inductiva.client.apis.paths.version_check import VersionCheck
from inductiva.client.apis.paths.users_quotas import UsersQuotas
from inductiva.client.apis.paths.users_info import UsersInfo
from inductiva.client.apis.paths.users_capabilities import UsersCapabilities
from inductiva.client.apis.paths.users_costs import UsersCosts
from inductiva.client.apis.paths.users_organization_costs import UsersOrganizationCosts
from inductiva.client.apis.paths.projects import Projects
from inductiva.client.apis.paths.projects_name import ProjectsName
from inductiva.client.apis.paths.metrics_users_username_activity import MetricsUsersUsernameActivity
from inductiva.client.apis.paths.metrics_users_username_cost_over_time import MetricsUsersUsernameCostOverTime
from inductiva.client.apis.paths.metrics_users_username_task_status_overview import MetricsUsersUsernameTaskStatusOverview
from inductiva.client.apis.paths.metrics_users_username_computation_time_trend import MetricsUsersUsernameComputationTimeTrend
from inductiva.client.apis.paths.metrics_users_username_tasks_overview import MetricsUsersUsernameTasksOverview
from inductiva.client.apis.paths.metrics_users_username_most_used_machine_types import MetricsUsersUsernameMostUsedMachineTypes
from inductiva.client.apis.paths.metrics_users_username_most_used_simulators_overview import MetricsUsersUsernameMostUsedSimulatorsOverview
from inductiva.client.apis.paths.metrics_usage_statistics import MetricsUsageStatistics

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
        PathValues.TASKS_TASK_ID:
            TasksTaskId,
        PathValues.TASKS:
            Tasks,
        PathValues.TASKS_TASK_ID_STATUS:
            TasksTaskIdStatus,
        PathValues.TASKS_TASK_ID_OUTPUT_LIST:
            TasksTaskIdOutputList,
        PathValues.TASKS_TASK_ID_DOWNLOAD_INPUT_URL:
            TasksTaskIdDownloadInputUrl,
        PathValues.TASKS_TASK_ID_DOWNLOAD_OUTPUT_URL:
            TasksTaskIdDownloadOutputUrl,
        PathValues.TASKS_TASK_ID_RESUBMIT:
            TasksTaskIdResubmit,
        PathValues.TASKS_TASK_ID_KILL:
            TasksTaskIdKill,
        PathValues.TASKS_TASK_ID_DISABLE_LOGS:
            TasksTaskIdDisableLogs,
        PathValues.TASKS_TASK_ID_FILES:
            TasksTaskIdFiles,
        PathValues.TASKS_TASK_ID_REGISTER:
            TasksTaskIdRegister,
        PathValues.TASKS_TASK_ID_OFFER:
            TasksTaskIdOffer,
        PathValues.TASKS_TASK_ID_MESSAGE:
            TasksTaskIdMessage,
        PathValues.ADMIN_USERS:
            AdminUsers,
        PathValues.ADMIN_USERS_EMAIL_TERMS_AND_CONDITIONS:
            AdminUsersEmailTermsAndConditions,
        PathValues.ADMIN_USERS_USERNAME_ORGANIZATION:
            AdminUsersUsernameOrganization,
        PathValues.ADMIN_USERS_USERNAME_TIER:
            AdminUsersUsernameTier,
        PathValues.ADMIN_USERS_USERNAME_CREDITS:
            AdminUsersUsernameCredits,
        PathValues.ADMIN_USERS_EMAIL_API_KEY:
            AdminUsersEmailApiKey,
        PathValues.ADMIN_USERS_EMAIL:
            AdminUsersEmail,
        PathValues.ADMIN_USERS_USERNAME_STORAGE_SIZE_FS:
            AdminUsersUsernameStorageSizeFs,
        PathValues.ADMIN_USERS_USERNAME_STORAGE_SIZE:
            AdminUsersUsernameStorageSize,
        PathValues.ADMIN_USERS_USERNAME_TASKS:
            AdminUsersUsernameTasks,
        PathValues.ADMIN_USERS_USERNAME_CAPABILITIES:
            AdminUsersUsernameCapabilities,
        PathValues.ADMIN_GROUPS:
            AdminGroups,
        PathValues.ADMIN_GROUPS_ACTIVE:
            AdminGroupsActive,
        PathValues.ADMIN_ACTIVE_TASKS:
            AdminActiveTasks,
        PathValues.ADMIN_GROUPS_MACHINE_GROUP_ID_TERMINATE:
            AdminGroupsMachineGroupIdTerminate,
        PathValues.ADMIN_ORGANIZATIONS:
            AdminOrganizations,
        PathValues.ADMIN_ORGANIZATIONS_ORGANIZATION_ID:
            AdminOrganizationsOrganizationId,
        PathValues.ADMIN_ORGANIZATIONS_COSTS:
            AdminOrganizationsCosts,
        PathValues.ADMIN_TIERS:
            AdminTiers,
        PathValues.ADMIN_TERMINATE_MACHINE_GROUPS_CREDITS_EXHAUSTED:
            AdminTerminateMachineGroupsCreditsExhausted,
        PathValues.TASKRUNNER_REGISTER:
            TaskRunnerRegister,
        PathValues.TASKRUNNER_MACHINE_ID:
            TaskRunnerMachineId,
        PathValues.TASKRUNNER_MACHINE_ID_TASK:
            TaskRunnerMachineIdTask,
        PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_MESSAGE:
            TaskRunnerMachineIdTaskTaskIdMessage,
        PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_MESSAGE_UNBLOCK:
            TaskRunnerMachineIdTaskTaskIdMessageUnblock,
        PathValues.TASKRUNNER_MACHINE_ID_EVENT:
            TaskRunnerMachineIdEvent,
        PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_OPERATION:
            TaskRunnerMachineIdTaskTaskIdOperation,
        PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_OPERATION_OPERATION_ID_DONE:
            TaskRunnerMachineIdTaskTaskIdOperationOperationIdDone,
        PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_DOWNLOAD_INPUT_URL:
            TaskRunnerMachineIdTaskTaskIdDownloadInputUrl,
        PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_UPLOAD_OUTPUT_URL:
            TaskRunnerMachineIdTaskTaskIdUploadOutputUrl,
        PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_METRIC:
            TaskRunnerMachineIdTaskTaskIdMetric,
        PathValues.TASKRUNNER_MACHINE_ID_RESIZE_DISK:
            TaskRunnerMachineIdResizeDisk,
        PathValues.TASKRUNNER_MACHINE_ID_RESIZE_DISK_DONE:
            TaskRunnerMachineIdResizeDiskDone,
        PathValues.TASKRUNNER_MACHINE_ID_DOWNLOAD_URLS:
            TaskRunnerMachineIdDownloadUrls,
        PathValues.COMPUTE_GROUP:
            ComputeGroup,
        PathValues.COMPUTE_GROUP_START:
            ComputeGroupStart,
        PathValues.COMPUTE_PRICE:
            ComputePrice,
        PathValues.COMPUTE_GROUPS:
            ComputeGroups,
        PathValues.COMPUTE_GROUPS_HISTORY:
            ComputeGroupsHistory,
        PathValues.COMPUTE_GROUP_STATUS:
            ComputeGroupStatus,
        PathValues.COMPUTE_MACHINE_TYPES:
            ComputeMachineTypes,
        PathValues.COMPUTE_GROUP_NAME:
            ComputeGroupName,
        PathValues.STORAGE_SIZE:
            StorageSize,
        PathValues.STORAGE_COST:
            StorageCost,
        PathValues.STORAGE_CONTENTS:
            StorageContents,
        PathValues.STORAGE_INPUT_URL:
            StorageInputUrl,
        PathValues.STORAGE_SIGNEDURLS:
            StorageSignedUrls,
        PathValues.STORAGE_INPUT_NOTIFY:
            StorageInputNotify,
        PathValues.STORAGE_INPUT_REMOTE:
            StorageInputRemote,
        PathValues.STORAGE_:
            Storage,
        PathValues.STORAGE_EXPORT:
            StorageExport,
        PathValues.STORAGE_OPERATIONS_OPERATION_ID:
            StorageOperationsOperationId,
        PathValues.STORAGE_OPERATIONS:
            StorageOperations,
        PathValues.VERSION:
            Version,
        PathValues.VERSIONCHECK:
            VersionCheck,
        PathValues.USERS_QUOTAS:
            UsersQuotas,
        PathValues.USERS_INFO:
            UsersInfo,
        PathValues.USERS_CAPABILITIES:
            UsersCapabilities,
        PathValues.USERS_COSTS:
            UsersCosts,
        PathValues.USERS_ORGANIZATION_COSTS:
            UsersOrganizationCosts,
        PathValues.PROJECTS:
            Projects,
        PathValues.PROJECTS_NAME:
            ProjectsName,
        PathValues.METRICS_USERS_USERNAME_ACTIVITY:
            MetricsUsersUsernameActivity,
        PathValues.METRICS_USERS_USERNAME_COST_OVER_TIME:
            MetricsUsersUsernameCostOverTime,
        PathValues.METRICS_USERS_USERNAME_TASK_STATUS_OVERVIEW:
            MetricsUsersUsernameTaskStatusOverview,
        PathValues.METRICS_USERS_USERNAME_COMPUTATION_TIME_TREND:
            MetricsUsersUsernameComputationTimeTrend,
        PathValues.METRICS_USERS_USERNAME_TASKS_OVERVIEW:
            MetricsUsersUsernameTasksOverview,
        PathValues.METRICS_USERS_USERNAME_MOST_USED_MACHINE_TYPES:
            MetricsUsersUsernameMostUsedMachineTypes,
        PathValues.METRICS_USERS_USERNAME_MOST_USED_SIMULATORS_OVERVIEW:
            MetricsUsersUsernameMostUsedSimulatorsOverview,
        PathValues.METRICS_USAGE_STATISTICS:
            MetricsUsageStatistics,
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
    PathValues.TASKS_TASK_ID:
        TasksTaskId,
    PathValues.TASKS:
        Tasks,
    PathValues.TASKS_TASK_ID_STATUS:
        TasksTaskIdStatus,
    PathValues.TASKS_TASK_ID_OUTPUT_LIST:
        TasksTaskIdOutputList,
    PathValues.TASKS_TASK_ID_DOWNLOAD_INPUT_URL:
        TasksTaskIdDownloadInputUrl,
    PathValues.TASKS_TASK_ID_DOWNLOAD_OUTPUT_URL:
        TasksTaskIdDownloadOutputUrl,
    PathValues.TASKS_TASK_ID_RESUBMIT:
        TasksTaskIdResubmit,
    PathValues.TASKS_TASK_ID_KILL:
        TasksTaskIdKill,
    PathValues.TASKS_TASK_ID_DISABLE_LOGS:
        TasksTaskIdDisableLogs,
    PathValues.TASKS_TASK_ID_FILES:
        TasksTaskIdFiles,
    PathValues.TASKS_TASK_ID_REGISTER:
        TasksTaskIdRegister,
    PathValues.TASKS_TASK_ID_OFFER:
        TasksTaskIdOffer,
    PathValues.TASKS_TASK_ID_MESSAGE:
        TasksTaskIdMessage,
    PathValues.ADMIN_USERS:
        AdminUsers,
    PathValues.ADMIN_USERS_EMAIL_TERMS_AND_CONDITIONS:
        AdminUsersEmailTermsAndConditions,
    PathValues.ADMIN_USERS_USERNAME_ORGANIZATION:
        AdminUsersUsernameOrganization,
    PathValues.ADMIN_USERS_USERNAME_TIER:
        AdminUsersUsernameTier,
    PathValues.ADMIN_USERS_USERNAME_CREDITS:
        AdminUsersUsernameCredits,
    PathValues.ADMIN_USERS_EMAIL_API_KEY:
        AdminUsersEmailApiKey,
    PathValues.ADMIN_USERS_EMAIL:
        AdminUsersEmail,
    PathValues.ADMIN_USERS_USERNAME_STORAGE_SIZE_FS:
        AdminUsersUsernameStorageSizeFs,
    PathValues.ADMIN_USERS_USERNAME_STORAGE_SIZE:
        AdminUsersUsernameStorageSize,
    PathValues.ADMIN_USERS_USERNAME_TASKS:
        AdminUsersUsernameTasks,
    PathValues.ADMIN_USERS_USERNAME_CAPABILITIES:
        AdminUsersUsernameCapabilities,
    PathValues.ADMIN_GROUPS:
        AdminGroups,
    PathValues.ADMIN_GROUPS_ACTIVE:
        AdminGroupsActive,
    PathValues.ADMIN_ACTIVE_TASKS:
        AdminActiveTasks,
    PathValues.ADMIN_GROUPS_MACHINE_GROUP_ID_TERMINATE:
        AdminGroupsMachineGroupIdTerminate,
    PathValues.ADMIN_ORGANIZATIONS:
        AdminOrganizations,
    PathValues.ADMIN_ORGANIZATIONS_ORGANIZATION_ID:
        AdminOrganizationsOrganizationId,
    PathValues.ADMIN_ORGANIZATIONS_COSTS:
        AdminOrganizationsCosts,
    PathValues.ADMIN_TIERS:
        AdminTiers,
    PathValues.ADMIN_TERMINATE_MACHINE_GROUPS_CREDITS_EXHAUSTED:
        AdminTerminateMachineGroupsCreditsExhausted,
    PathValues.TASKRUNNER_REGISTER:
        TaskRunnerRegister,
    PathValues.TASKRUNNER_MACHINE_ID:
        TaskRunnerMachineId,
    PathValues.TASKRUNNER_MACHINE_ID_TASK:
        TaskRunnerMachineIdTask,
    PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_MESSAGE:
        TaskRunnerMachineIdTaskTaskIdMessage,
    PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_MESSAGE_UNBLOCK:
        TaskRunnerMachineIdTaskTaskIdMessageUnblock,
    PathValues.TASKRUNNER_MACHINE_ID_EVENT:
        TaskRunnerMachineIdEvent,
    PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_OPERATION:
        TaskRunnerMachineIdTaskTaskIdOperation,
    PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_OPERATION_OPERATION_ID_DONE:
        TaskRunnerMachineIdTaskTaskIdOperationOperationIdDone,
    PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_DOWNLOAD_INPUT_URL:
        TaskRunnerMachineIdTaskTaskIdDownloadInputUrl,
    PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_UPLOAD_OUTPUT_URL:
        TaskRunnerMachineIdTaskTaskIdUploadOutputUrl,
    PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_METRIC:
        TaskRunnerMachineIdTaskTaskIdMetric,
    PathValues.TASKRUNNER_MACHINE_ID_RESIZE_DISK:
        TaskRunnerMachineIdResizeDisk,
    PathValues.TASKRUNNER_MACHINE_ID_RESIZE_DISK_DONE:
        TaskRunnerMachineIdResizeDiskDone,
    PathValues.TASKRUNNER_MACHINE_ID_DOWNLOAD_URLS:
        TaskRunnerMachineIdDownloadUrls,
    PathValues.COMPUTE_GROUP:
        ComputeGroup,
    PathValues.COMPUTE_GROUP_START:
        ComputeGroupStart,
    PathValues.COMPUTE_PRICE:
        ComputePrice,
    PathValues.COMPUTE_GROUPS:
        ComputeGroups,
    PathValues.COMPUTE_GROUPS_HISTORY:
        ComputeGroupsHistory,
    PathValues.COMPUTE_GROUP_STATUS:
        ComputeGroupStatus,
    PathValues.COMPUTE_MACHINE_TYPES:
        ComputeMachineTypes,
    PathValues.COMPUTE_GROUP_NAME:
        ComputeGroupName,
    PathValues.STORAGE_SIZE:
        StorageSize,
    PathValues.STORAGE_COST:
        StorageCost,
    PathValues.STORAGE_CONTENTS:
        StorageContents,
    PathValues.STORAGE_INPUT_URL:
        StorageInputUrl,
    PathValues.STORAGE_SIGNEDURLS:
        StorageSignedUrls,
    PathValues.STORAGE_INPUT_NOTIFY:
        StorageInputNotify,
    PathValues.STORAGE_INPUT_REMOTE:
        StorageInputRemote,
    PathValues.STORAGE_:
        Storage,
    PathValues.STORAGE_EXPORT:
        StorageExport,
    PathValues.STORAGE_OPERATIONS_OPERATION_ID:
        StorageOperationsOperationId,
    PathValues.STORAGE_OPERATIONS:
        StorageOperations,
    PathValues.VERSION:
        Version,
    PathValues.VERSIONCHECK:
        VersionCheck,
    PathValues.USERS_QUOTAS:
        UsersQuotas,
    PathValues.USERS_INFO:
        UsersInfo,
    PathValues.USERS_CAPABILITIES:
        UsersCapabilities,
    PathValues.USERS_COSTS:
        UsersCosts,
    PathValues.USERS_ORGANIZATION_COSTS:
        UsersOrganizationCosts,
    PathValues.PROJECTS:
        Projects,
    PathValues.PROJECTS_NAME:
        ProjectsName,
    PathValues.METRICS_USERS_USERNAME_ACTIVITY:
        MetricsUsersUsernameActivity,
    PathValues.METRICS_USERS_USERNAME_COST_OVER_TIME:
        MetricsUsersUsernameCostOverTime,
    PathValues.METRICS_USERS_USERNAME_TASK_STATUS_OVERVIEW:
        MetricsUsersUsernameTaskStatusOverview,
    PathValues.METRICS_USERS_USERNAME_COMPUTATION_TIME_TREND:
        MetricsUsersUsernameComputationTimeTrend,
    PathValues.METRICS_USERS_USERNAME_TASKS_OVERVIEW:
        MetricsUsersUsernameTasksOverview,
    PathValues.METRICS_USERS_USERNAME_MOST_USED_MACHINE_TYPES:
        MetricsUsersUsernameMostUsedMachineTypes,
    PathValues.METRICS_USERS_USERNAME_MOST_USED_SIMULATORS_OVERVIEW:
        MetricsUsersUsernameMostUsedSimulatorsOverview,
    PathValues.METRICS_USAGE_STATISTICS:
        MetricsUsageStatistics,
})
