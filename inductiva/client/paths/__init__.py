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
    TASKS_TASK_ID_STDOUT_TAIL = "/tasks/{task_id}/stdout_tail"
    TASKS_TASK_ID_RESOURCES_TAIL = "/tasks/{task_id}/resources_tail"
    ADMIN_USERS = "/admin/users"
    ADMIN_USERS_EMAIL_API_KEY = "/admin/users/{email}/api_key"
    ADMIN_USERS_USERNAME_TASKS = "/admin/users/{username}/tasks"
    ADMIN_GROUPS = "/admin/groups"
    EXECUTERTRACKER_REGISTER = "/executer-tracker/register"
    COMPUTE_GROUP = "/compute/group"
    COMPUTE_GROUP_START = "/compute/group/start"
    COMPUTE_GROUP_ELASTIC = "/compute/group/elastic"
    COMPUTE_GROUP_GPU = "/compute/group/gpu"
    COMPUTE_PRICE = "/compute/price"
    COMPUTE_STATUS = "/compute/status"
    COMPUTE_GROUP_STATUS = "/compute/group_status"
    COMPUTE_GROUPS = "/compute/groups"
    COMPUTE_GROUP_NAME = "/compute/group/{name}"
    STORAGE_SIZE = "/storage/size"
    STORAGE_CONTENTS = "/storage/contents"
    STORAGE_DIR_NAME = "/storage/{dir_name}"
    VERSION = "/version"
    VERSIONCHECK = "/version-check"
