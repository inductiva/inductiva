import typing_extensions

from inductiva.client.paths import PathValues
from inductiva.client.apis.paths.tasks_submit import TasksSubmit
from inductiva.client.apis.paths.tasks_task_id_input_uploaded import TasksTaskIdInputUploaded
from inductiva.client.apis.paths.tasks_task_id import TasksTaskId
from inductiva.client.apis.paths.tasks import Tasks
from inductiva.client.apis.paths.tasks_task_id_status import TasksTaskIdStatus
from inductiva.client.apis.paths.tasks_task_id_kill import TasksTaskIdKill
from inductiva.client.apis.paths.tasks_task_id_register import TasksTaskIdRegister
from inductiva.client.apis.paths.tasks_task_id_offer import TasksTaskIdOffer
from inductiva.client.apis.paths.tasks_task_id_message import TasksTaskIdMessage
from inductiva.client.apis.paths.tasks_task_id_metadata import TasksTaskIdMetadata
from inductiva.client.apis.paths.admin_users import AdminUsers
from inductiva.client.apis.paths.admin_users_email_terms_and_conditions import AdminUsersEmailTermsAndConditions
from inductiva.client.apis.paths.admin_users_username_organization import AdminUsersUsernameOrganization
from inductiva.client.apis.paths.admin_users_username_credits import AdminUsersUsernameCredits
from inductiva.client.apis.paths.admin_user_emails import AdminUserEmails
from inductiva.client.apis.paths.admin_user_buckets import AdminUserBuckets
from inductiva.client.apis.paths.admin_users_email_api_key import AdminUsersEmailApiKey
from inductiva.client.apis.paths.admin_users_email import AdminUsersEmail
from inductiva.client.apis.paths.admin_users_username_storage_size_fs import AdminUsersUsernameStorageSizeFs
from inductiva.client.apis.paths.admin_users_username_storage_size import AdminUsersUsernameStorageSize
from inductiva.client.apis.paths.admin_users_username_tasks import AdminUsersUsernameTasks
from inductiva.client.apis.paths.admin_users_username_capabilities import AdminUsersUsernameCapabilities
from inductiva.client.apis.paths.admin_groups import AdminGroups
from inductiva.client.apis.paths.admin_groups_active import AdminGroupsActive
from inductiva.client.apis.paths.admin_groups_id import AdminGroupsId
from inductiva.client.apis.paths.admin_active_tasks import AdminActiveTasks
from inductiva.client.apis.paths.admin_groups_machine_group_id_terminate import AdminGroupsMachineGroupIdTerminate
from inductiva.client.apis.paths.admin_organizations import AdminOrganizations
from inductiva.client.apis.paths.admin_organizations_organization_id import AdminOrganizationsOrganizationId
from inductiva.client.apis.paths.admin_organizations_costs import AdminOrganizationsCosts
from inductiva.client.apis.paths.admin_tiers import AdminTiers
from inductiva.client.apis.paths.admin_terminate_machine_groups_credits_exhausted import AdminTerminateMachineGroupsCreditsExhausted
from inductiva.client.apis.paths.admin_import_provider_costs import AdminImportProviderCosts
from inductiva.client.apis.paths.admin_machine_id_event import AdminMachineIdEvent
from inductiva.client.apis.paths.admin_users_username_costs_fee_percentage import AdminUsersUsernameCostsFeePercentage
from inductiva.client.apis.paths.admin_organizations_organization_id_costs_fee_percentage import AdminOrganizationsOrganizationIdCostsFeePercentage
from inductiva.client.apis.paths.admin_organizations_organization_id_terminate_resources_credits_threshold import AdminOrganizationsOrganizationIdTerminateResourcesCreditsThreshold
from inductiva.client.apis.paths.admin_top_ups import AdminTopUps
from inductiva.client.apis.paths.admin_tasks import AdminTasks
from inductiva.client.apis.paths.admin_users_email_stripe_customer_id import AdminUsersEmailStripeCustomerId
from inductiva.client.apis.paths.admin_feature_flags_name import AdminFeatureFlagsName
from inductiva.client.apis.paths.admin_feature_flags_ import AdminFeatureFlags
from inductiva.client.apis.paths.admin_alerts_check_credits import AdminAlertsCheckCredits
from inductiva.client.apis.paths.admin_alerts_check_tasks import AdminAlertsCheckTasks
from inductiva.client.apis.paths.task_runner_register import TaskRunnerRegister
from inductiva.client.apis.paths.task_runner_machine_id import TaskRunnerMachineId
from inductiva.client.apis.paths.task_runner_machine_id_task import TaskRunnerMachineIdTask
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_message import TaskRunnerMachineIdTaskTaskIdMessage
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_message_unblock import TaskRunnerMachineIdTaskTaskIdMessageUnblock
from inductiva.client.apis.paths.task_runner_machine_id_event import TaskRunnerMachineIdEvent
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_operation import TaskRunnerMachineIdTaskTaskIdOperation
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_operation_operation_id_done import TaskRunnerMachineIdTaskTaskIdOperationOperationIdDone
from inductiva.client.apis.paths.task_runner_machine_id_task_task_id_metric import TaskRunnerMachineIdTaskTaskIdMetric
from inductiva.client.apis.paths.task_runner_machine_id_resize_disk import TaskRunnerMachineIdResizeDisk
from inductiva.client.apis.paths.task_runner_machine_id_resize_disk_done import TaskRunnerMachineIdResizeDiskDone
from inductiva.client.apis.paths.compute_group import ComputeGroup
from inductiva.client.apis.paths.compute_group_start import ComputeGroupStart
from inductiva.client.apis.paths.compute_price import ComputePrice
from inductiva.client.apis.paths.compute_groups import ComputeGroups
from inductiva.client.apis.paths.compute_groups_history import ComputeGroupsHistory
from inductiva.client.apis.paths.compute_machine_types import ComputeMachineTypes
from inductiva.client.apis.paths.compute_group_name import ComputeGroupName
from inductiva.client.apis.paths.compute_group_machine_group_id_sharing import ComputeGroupMachineGroupIdSharing
from inductiva.client.apis.paths.storage_size import StorageSize
from inductiva.client.apis.paths.storage_cost import StorageCost
from inductiva.client.apis.paths.storage_contents import StorageContents
from inductiva.client.apis.paths.storage_zip_contents import StorageZipContents
from inductiva.client.apis.paths.storage_signed_urls import StorageSignedUrls
from inductiva.client.apis.paths.storage_input_remote import StorageInputRemote
from inductiva.client.apis.paths.storage_ import Storage
from inductiva.client.apis.paths.storage_copy import StorageCopy
from inductiva.client.apis.paths.storage_operations_operation_id import StorageOperationsOperationId
from inductiva.client.apis.paths.storage_operations import StorageOperations
from inductiva.client.apis.paths.storage_export_multipart import StorageExportMultipart
from inductiva.client.apis.paths.storage_update_operation_status import StorageUpdateOperationStatus
from inductiva.client.apis.paths.version import Version
from inductiva.client.apis.paths.version_check import VersionCheck
from inductiva.client.apis.paths.users_quotas import UsersQuotas
from inductiva.client.apis.paths.users_info import UsersInfo
from inductiva.client.apis.paths.users_capabilities import UsersCapabilities
from inductiva.client.apis.paths.users_costs import UsersCosts
from inductiva.client.apis.paths.users_organization_costs import UsersOrganizationCosts
from inductiva.client.apis.paths.users_top_ups import UsersTopUps
from inductiva.client.apis.paths.projects import Projects
from inductiva.client.apis.paths.projects_name import ProjectsName
from inductiva.client.apis.paths.projects_name_task_task_id_add import ProjectsNameTaskTaskIdAdd
from inductiva.client.apis.paths.projects_name_task_task_id_remove import ProjectsNameTaskTaskIdRemove
from inductiva.client.apis.paths.projects_name_metadata import ProjectsNameMetadata
from inductiva.client.apis.paths.pubsub_notify_file_change import PubsubNotifyFileChange
from inductiva.client.apis.paths.metrics_users_username_activity import MetricsUsersUsernameActivity
from inductiva.client.apis.paths.metrics_users_username_cost_over_time import MetricsUsersUsernameCostOverTime
from inductiva.client.apis.paths.metrics_users_username_task_status_overview import MetricsUsersUsernameTaskStatusOverview
from inductiva.client.apis.paths.metrics_users_username_computation_time_trend import MetricsUsersUsernameComputationTimeTrend
from inductiva.client.apis.paths.metrics_users_username_tasks_overview import MetricsUsersUsernameTasksOverview
from inductiva.client.apis.paths.metrics_users_username_most_used_machine_types import MetricsUsersUsernameMostUsedMachineTypes
from inductiva.client.apis.paths.metrics_users_username_most_used_simulators_overview import MetricsUsersUsernameMostUsedSimulatorsOverview
from inductiva.client.apis.paths.metrics_usage_statistics import MetricsUsageStatistics
from inductiva.client.apis.paths.events_ import Events
from inductiva.client.apis.paths.events_event_id import EventsEventId
from inductiva.client.apis.paths.simulators_available_images import SimulatorsAvailableImages

PathToApi = typing_extensions.TypedDict(
    'PathToApi', {
        PathValues.TASKS_SUBMIT:
            TasksSubmit,
        PathValues.TASKS_TASK_ID_INPUT_UPLOADED:
            TasksTaskIdInputUploaded,
        PathValues.TASKS_TASK_ID:
            TasksTaskId,
        PathValues.TASKS:
            Tasks,
        PathValues.TASKS_TASK_ID_STATUS:
            TasksTaskIdStatus,
        PathValues.TASKS_TASK_ID_KILL:
            TasksTaskIdKill,
        PathValues.TASKS_TASK_ID_REGISTER:
            TasksTaskIdRegister,
        PathValues.TASKS_TASK_ID_OFFER:
            TasksTaskIdOffer,
        PathValues.TASKS_TASK_ID_MESSAGE:
            TasksTaskIdMessage,
        PathValues.TASKS_TASK_ID_METADATA:
            TasksTaskIdMetadata,
        PathValues.ADMIN_USERS:
            AdminUsers,
        PathValues.ADMIN_USERS_EMAIL_TERMS_AND_CONDITIONS:
            AdminUsersEmailTermsAndConditions,
        PathValues.ADMIN_USERS_USERNAME_ORGANIZATION:
            AdminUsersUsernameOrganization,
        PathValues.ADMIN_USERS_USERNAME_CREDITS:
            AdminUsersUsernameCredits,
        PathValues.ADMIN_USER_EMAILS:
            AdminUserEmails,
        PathValues.ADMIN_USER_BUCKETS:
            AdminUserBuckets,
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
        PathValues.ADMIN_GROUPS_ID:
            AdminGroupsId,
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
        PathValues.ADMIN_IMPORT_PROVIDER_COSTS:
            AdminImportProviderCosts,
        PathValues.ADMIN_MACHINE_ID_EVENT:
            AdminMachineIdEvent,
        PathValues.ADMIN_USERS_USERNAME_COSTS_FEE_PERCENTAGE:
            AdminUsersUsernameCostsFeePercentage,
        PathValues.ADMIN_ORGANIZATIONS_ORGANIZATION_ID_COSTS_FEE_PERCENTAGE:
            AdminOrganizationsOrganizationIdCostsFeePercentage,
        PathValues.ADMIN_ORGANIZATIONS_ORGANIZATION_ID_TERMINATE_RESOURCES_CREDITS_THRESHOLD:
            AdminOrganizationsOrganizationIdTerminateResourcesCreditsThreshold,
        PathValues.ADMIN_TOPUPS:
            AdminTopUps,
        PathValues.ADMIN_TASKS:
            AdminTasks,
        PathValues.ADMIN_USERS_EMAIL_STRIPE_CUSTOMER_ID:
            AdminUsersEmailStripeCustomerId,
        PathValues.ADMIN_FEATUREFLAGS_NAME:
            AdminFeatureFlagsName,
        PathValues.ADMIN_FEATUREFLAGS_:
            AdminFeatureFlags,
        PathValues.ADMIN_ALERTS_CHECK_CREDITS:
            AdminAlertsCheckCredits,
        PathValues.ADMIN_ALERTS_CHECK_TASKS:
            AdminAlertsCheckTasks,
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
        PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_METRIC:
            TaskRunnerMachineIdTaskTaskIdMetric,
        PathValues.TASKRUNNER_MACHINE_ID_RESIZE_DISK:
            TaskRunnerMachineIdResizeDisk,
        PathValues.TASKRUNNER_MACHINE_ID_RESIZE_DISK_DONE:
            TaskRunnerMachineIdResizeDiskDone,
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
        PathValues.COMPUTE_MACHINE_TYPES:
            ComputeMachineTypes,
        PathValues.COMPUTE_GROUP_NAME:
            ComputeGroupName,
        PathValues.COMPUTE_GROUP_MACHINE_GROUP_ID_SHARING:
            ComputeGroupMachineGroupIdSharing,
        PathValues.STORAGE_SIZE:
            StorageSize,
        PathValues.STORAGE_COST:
            StorageCost,
        PathValues.STORAGE_CONTENTS:
            StorageContents,
        PathValues.STORAGE_ZIPCONTENTS:
            StorageZipContents,
        PathValues.STORAGE_SIGNEDURLS:
            StorageSignedUrls,
        PathValues.STORAGE_INPUT_REMOTE:
            StorageInputRemote,
        PathValues.STORAGE_:
            Storage,
        PathValues.STORAGE_COPY:
            StorageCopy,
        PathValues.STORAGE_OPERATIONS_OPERATION_ID:
            StorageOperationsOperationId,
        PathValues.STORAGE_OPERATIONS:
            StorageOperations,
        PathValues.STORAGE_EXPORT_MULTIPART:
            StorageExportMultipart,
        PathValues.STORAGE_UPDATE_OPERATION_STATUS:
            StorageUpdateOperationStatus,
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
        PathValues.USERS_TOPUPS:
            UsersTopUps,
        PathValues.PROJECTS:
            Projects,
        PathValues.PROJECTS_NAME:
            ProjectsName,
        PathValues.PROJECTS_NAME_TASK_TASK_ID_ADD:
            ProjectsNameTaskTaskIdAdd,
        PathValues.PROJECTS_NAME_TASK_TASK_ID_REMOVE:
            ProjectsNameTaskTaskIdRemove,
        PathValues.PROJECTS_NAME_METADATA:
            ProjectsNameMetadata,
        PathValues.PUBSUB_NOTIFY_FILE_CHANGE:
            PubsubNotifyFileChange,
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
        PathValues.EVENTS_:
            Events,
        PathValues.EVENTS_EVENT_ID:
            EventsEventId,
        PathValues.SIMULATORS_AVAILABLEIMAGES:
            SimulatorsAvailableImages,
    })

path_to_api = PathToApi({
    PathValues.TASKS_SUBMIT:
        TasksSubmit,
    PathValues.TASKS_TASK_ID_INPUT_UPLOADED:
        TasksTaskIdInputUploaded,
    PathValues.TASKS_TASK_ID:
        TasksTaskId,
    PathValues.TASKS:
        Tasks,
    PathValues.TASKS_TASK_ID_STATUS:
        TasksTaskIdStatus,
    PathValues.TASKS_TASK_ID_KILL:
        TasksTaskIdKill,
    PathValues.TASKS_TASK_ID_REGISTER:
        TasksTaskIdRegister,
    PathValues.TASKS_TASK_ID_OFFER:
        TasksTaskIdOffer,
    PathValues.TASKS_TASK_ID_MESSAGE:
        TasksTaskIdMessage,
    PathValues.TASKS_TASK_ID_METADATA:
        TasksTaskIdMetadata,
    PathValues.ADMIN_USERS:
        AdminUsers,
    PathValues.ADMIN_USERS_EMAIL_TERMS_AND_CONDITIONS:
        AdminUsersEmailTermsAndConditions,
    PathValues.ADMIN_USERS_USERNAME_ORGANIZATION:
        AdminUsersUsernameOrganization,
    PathValues.ADMIN_USERS_USERNAME_CREDITS:
        AdminUsersUsernameCredits,
    PathValues.ADMIN_USER_EMAILS:
        AdminUserEmails,
    PathValues.ADMIN_USER_BUCKETS:
        AdminUserBuckets,
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
    PathValues.ADMIN_GROUPS_ID:
        AdminGroupsId,
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
    PathValues.ADMIN_IMPORT_PROVIDER_COSTS:
        AdminImportProviderCosts,
    PathValues.ADMIN_MACHINE_ID_EVENT:
        AdminMachineIdEvent,
    PathValues.ADMIN_USERS_USERNAME_COSTS_FEE_PERCENTAGE:
        AdminUsersUsernameCostsFeePercentage,
    PathValues.ADMIN_ORGANIZATIONS_ORGANIZATION_ID_COSTS_FEE_PERCENTAGE:
        AdminOrganizationsOrganizationIdCostsFeePercentage,
    PathValues.ADMIN_ORGANIZATIONS_ORGANIZATION_ID_TERMINATE_RESOURCES_CREDITS_THRESHOLD:
        AdminOrganizationsOrganizationIdTerminateResourcesCreditsThreshold,
    PathValues.ADMIN_TOPUPS:
        AdminTopUps,
    PathValues.ADMIN_TASKS:
        AdminTasks,
    PathValues.ADMIN_USERS_EMAIL_STRIPE_CUSTOMER_ID:
        AdminUsersEmailStripeCustomerId,
    PathValues.ADMIN_FEATUREFLAGS_NAME:
        AdminFeatureFlagsName,
    PathValues.ADMIN_FEATUREFLAGS_:
        AdminFeatureFlags,
    PathValues.ADMIN_ALERTS_CHECK_CREDITS:
        AdminAlertsCheckCredits,
    PathValues.ADMIN_ALERTS_CHECK_TASKS:
        AdminAlertsCheckTasks,
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
    PathValues.TASKRUNNER_MACHINE_ID_TASK_TASK_ID_METRIC:
        TaskRunnerMachineIdTaskTaskIdMetric,
    PathValues.TASKRUNNER_MACHINE_ID_RESIZE_DISK:
        TaskRunnerMachineIdResizeDisk,
    PathValues.TASKRUNNER_MACHINE_ID_RESIZE_DISK_DONE:
        TaskRunnerMachineIdResizeDiskDone,
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
    PathValues.COMPUTE_MACHINE_TYPES:
        ComputeMachineTypes,
    PathValues.COMPUTE_GROUP_NAME:
        ComputeGroupName,
    PathValues.COMPUTE_GROUP_MACHINE_GROUP_ID_SHARING:
        ComputeGroupMachineGroupIdSharing,
    PathValues.STORAGE_SIZE:
        StorageSize,
    PathValues.STORAGE_COST:
        StorageCost,
    PathValues.STORAGE_CONTENTS:
        StorageContents,
    PathValues.STORAGE_ZIPCONTENTS:
        StorageZipContents,
    PathValues.STORAGE_SIGNEDURLS:
        StorageSignedUrls,
    PathValues.STORAGE_INPUT_REMOTE:
        StorageInputRemote,
    PathValues.STORAGE_:
        Storage,
    PathValues.STORAGE_COPY:
        StorageCopy,
    PathValues.STORAGE_OPERATIONS_OPERATION_ID:
        StorageOperationsOperationId,
    PathValues.STORAGE_OPERATIONS:
        StorageOperations,
    PathValues.STORAGE_EXPORT_MULTIPART:
        StorageExportMultipart,
    PathValues.STORAGE_UPDATE_OPERATION_STATUS:
        StorageUpdateOperationStatus,
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
    PathValues.USERS_TOPUPS:
        UsersTopUps,
    PathValues.PROJECTS:
        Projects,
    PathValues.PROJECTS_NAME:
        ProjectsName,
    PathValues.PROJECTS_NAME_TASK_TASK_ID_ADD:
        ProjectsNameTaskTaskIdAdd,
    PathValues.PROJECTS_NAME_TASK_TASK_ID_REMOVE:
        ProjectsNameTaskTaskIdRemove,
    PathValues.PROJECTS_NAME_METADATA:
        ProjectsNameMetadata,
    PathValues.PUBSUB_NOTIFY_FILE_CHANGE:
        PubsubNotifyFileChange,
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
    PathValues.EVENTS_:
        Events,
    PathValues.EVENTS_EVENT_ID:
        EventsEventId,
    PathValues.SIMULATORS_AVAILABLEIMAGES:
        SimulatorsAvailableImages,
})
