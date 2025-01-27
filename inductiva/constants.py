"""Constants that can be set for the Inductiva client."""
import os
import pathlib
import platform

LOGS_WEBSOCKET_URL = os.environ.get("INDUCTIVA_TASK_LOGS_URL",
                                    "wss://logs.inductiva.ai")

TASK_RUNNER_IMAGE = os.environ.get("INDUCTIVA_TASK_RUNNER_IMAGE",
                                   "inductiva/task-runner:main")

FILE_TRACKER_IMAGE = os.environ.get("INDUCTIVA_FILE_TRACKER_IMAGE",
                                    "inductiva/file-tracker:main")

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

system = platform.system()
if platform.system().lower() == "windows":
    HOME_DIR = pathlib.Path.home() / "AppData" / "inductiva"
elif platform.system().lower() in ["linux", "darwin"]:
    HOME_DIR = pathlib.Path.home() / ".inductiva"
else:
    raise RuntimeError(f"Current operating system {system} not supported.")

LOGS_FILE_PATH = HOME_DIR / "inductiva.log"
API_KEY_FILE_PATH = HOME_DIR / "api_key"

TASK_OUTPUT_ZIP = "output.zip"

LOGIN_MESSAGE = "Please login with `inductiva auth login`."
