# do not import all endpoints into this module because that uses a lot of memory and stack frames
# if you need the ability to import all endpoints from this module, import them with
# from client.paths.task_task_id_status import Api

from client.paths import PathValues

path = PathValues.TASK_TASK_ID_STATUS