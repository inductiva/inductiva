# do not import all endpoints into this module because that uses a lot of memory and stack frames
# if you need the ability to import all endpoints from this module, import them with
# from inductiva.client.apis.path_to_api import path_to_api

import enum


class PathValues(str, enum.Enum):
    TASK_SUBMIT = "/task/submit"
    TASK_TASK_ID_INPUT = "/task/{task_id}/input"
    TASK_TASK_ID_STATUS = "/task/{task_id}/status"
    TASK_TASK_ID_OUTPUT = "/task/{task_id}/output"
    TASK_TASK_ID_KILL = "/task/{task_id}/kill"
    ADMIN_USER = "/admin/user"
    ADMIN_USER_USERNAME_TASKS = "/admin/user/{username}/tasks"
