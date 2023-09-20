# do not import all endpoints into this module because that uses a lot of memory and stack frames
# if you need the ability to import all endpoints from this module, import them with
# from inductiva.client.apis.path_to_api import path_to_api

import enum


class PathValues(str, enum.Enum):
    TASKS_SUBMIT = "/tasks/submit"
    TASKS_TASK_ID_INPUT = "/tasks/{task_id}/input"
    TASKS_TASK_ID = "/tasks/{task_id}"
    TASKS = "/tasks"
    TASKS_TASK_ID_STATUS = "/tasks/{task_id}/status"
    TASKS_TASK_ID_OUTPUT_LIST = "/tasks/{task_id}/output/list"
    TASKS_TASK_ID_OUTPUT = "/tasks/{task_id}/output"
    TASKS_TASK_ID_KILL = "/tasks/{task_id}/kill"
    ADMIN_USERS = "/admin/users"
    ADMIN_USERS_USERNAME = "/admin/users/{username}"
    ADMIN_USERS_USERNAME_TASKS = "/admin/users/{username}/tasks"
    ADMIN_TASKS = "/admin/tasks"
    ADMIN_GROUPS = "/admin/groups"
    EXECUTERS_REGISTER = "/executers/register"
    GCP_INSTANCES_GROUP = "/gcp_instances/group"
    GCP_INSTANCES_GROUP_START = "/gcp_instances/group/start"
    GCP_INSTANCES_PRICE = "/gcp_instances/price"
    GCP_INSTANCES_STATUS = "/gcp_instances/status"
    GCP_INSTANCES_GROUP_STATUS = "/gcp_instances/group_status"
    GCP_INSTANCES_STORAGE = "/gcp_instances/storage"
    GCP_INSTANCES_GROUPS = "/gcp_instances/groups"
    GCP_INSTANCES_GROUP_NAME = "/gcp_instances/group/{name}"
