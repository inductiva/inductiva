"""Constants that can be set for the Inductiva client."""
import os
import pathlib

LOGS_WEBSOCKET_URL = os.environ.get("INDUCTIVA_TASK_LOGS_URL",
                                    "wss://logs.inductiva.ai")

DEFAULT_QUEUE_MACHINE_TYPE = "c2-standard-4"

TASK_KILL_MAX_API_REQUESTS = 5

TASK_KILL_RETRY_SLEEP_SEC = 1.0

MAX_CONFIRMATION_LINES = 5

LOADER_COMMAND_PREFIX = "cmd_"

LOADER_IGNORE_PREFIX = "_"

# when printing the stack trace, how many lines to show
EXCEPTIONS_MAX_TRACEBACK_DEPTH = 2

TASK_FAILED_LINES_TO_DUMP = 10

INDUCTIVA_LOGS_WAIT_SLEEP_TIME = 1.0

HOME_DIR = pathlib.Path.home() / ".inductiva"
