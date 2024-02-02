"""Constants that can be set for the Inductiva client."""
from inductiva.client import models
import os

LOGS_WEBSOCKET_URL = os.environ.get("INDUCTIVA_TASK_LOGS_URL",
                                    "wss://logs.inductiva.ai")

DEFAULT_QUEUE_MACHINE_TYPE = "c2-standard-4"

TASK_KILL_MAX_API_REQUESTS = 5

TASK_KILL_RETRY_SLEEP_SEC = 1.0

MAX_CONFIRMATION_LINES = 5

TASK_FAILED_STATUSES = {
    models.TaskStatusCode.FAILED, models.TaskStatusCode.KILLED,
    models.TaskStatusCode.EXECUTERFAILED,
    models.TaskStatusCode.EXECUTERTERMINATED,
    models.TaskStatusCode.EXECUTERTERMINATEDBYUSER,
    models.TaskStatusCode.SPOTINSTANCEPREEMPTED, models.TaskStatusCode.ZOMBIE
}

TASK_TERMINAL_STATUSES = {models.TaskStatusCode.SUCCESS
                         }.union(TASK_FAILED_STATUSES)

TASK_RUNNING_STATUSES = {
    models.TaskStatusCode.PENDINGINPUT, models.TaskStatusCode.SUBMITTED,
    models.TaskStatusCode.STARTED
}

LOADER_COMMAND_PREFIX = "cmd_"

LOADER_IGNORE_PREFIX = "_"
