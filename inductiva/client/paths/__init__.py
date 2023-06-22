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
    TASKS_TASK_ID_OUTPUT = "/tasks/{task_id}/output"
    TASKS_TASK_ID_KILL = "/tasks/{task_id}/kill"
    ADMIN_USERS = "/admin/users"
    ADMIN_USERS_USERNAME = "/admin/users/{username}"
    ADMIN_USERS_USERNAME_TASKS = "/admin/users/{username}/tasks"
    ADMIN_TASKS = "/admin/tasks"
    EXECUTERS_REGISTER = "/executers/register"
    EXECUTERS_POOLS = "/executers/pools"
    GCP_INSTANCES = "/gcp_instances"
    GCP_INSTANCES_GROUP = "/gcp_instances/group"
    GCP_INSTANCES_PRICE = "/gcp_instances/price"
