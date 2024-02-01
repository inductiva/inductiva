"""Constants that can be set for the Inductiva client."""
import os

LOGS_WEBSOCKET_URL = os.environ.get("INDUCTIVA_TASK_LOGS_URL",
                                    "wss://logs.inductiva.ai")

DEFAULT_QUEUE_MACHINE_TYPE = "c2-standard-4"

TASK_KILL_MAX_API_REQUESTS = 5

TASK_KILL_RETRY_SLEEP_SEC = 1.0

LOADER_LOAD_PREFIX = "cmd_"

LOADER_IGNORE_PREFIX = "_"
