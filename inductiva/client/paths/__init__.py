# do not import all endpoints into this module because that uses a lot of memory and stack frames
# if you need the ability to import all endpoints from this module, import them with
# from inductiva.client.apis.path_to_api import path_to_api

import enum


class PathValues(str, enum.Enum):
    TASKS_AUTH = "/tasks/auth"
    TASKS_SUBMIT = "/tasks/submit"
    TASKS_TASK_ID_INPUT = "/tasks/{task_id}/input"
    TASKS_TASK_ID = "/tasks/{task_id}"
    TASKS = "/tasks"
    TASKS_TASK_ID_STATUS = "/tasks/{task_id}/status"
    TASKS_TASK_ID_OUTPUT_LIST = "/tasks/{task_id}/output/list"
    TASKS_TASK_ID_OUTPUT = "/tasks/{task_id}/output"
    TASKS_TASK_ID_KILL = "/tasks/{task_id}/kill"
    ADMIN_USERS = "/admin/users"
    ADMIN_USERS_EMAIL_API_KEY = "/admin/users/{email}/api_key"
    ADMIN_USERS_USERNAME_TASKS = "/admin/users/{username}/tasks"
    ADMIN_GROUPS = "/admin/groups"
    ADMIN_ACTIVE_TASKS = "/admin/active_tasks"
    EXECUTERTRACKER_REGISTER = "/executer-tracker/register"
    COMPUTE_GROUP = "/compute/group"
    COMPUTE_TYPE = "/compute/type"
    COMPUTE_GROUP_START = "/compute/group/start"
    COMPUTE_PRICE = "/compute/price"
    COMPUTE_GROUPS = "/compute/groups"
    COMPUTE_GROUP_STATUS = "/compute/group_status"
    STORAGE_SIZE = "/storage/size"
    STORAGE_CONTENTS = "/storage/contents"
    VERSION = "/version"
    VERSIONCHECK = "/version-check"
