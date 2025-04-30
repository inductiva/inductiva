"""Constants that can be set for the Inductiva client."""
import os
import pathlib
import platform

TURN_SERVER_URL = os.environ.get("INDUCTIVA_TURN_SERVER_URL",
                                 "webrtc.inductiva.ai:3478")

TASK_RUNNER_IMAGE = os.environ.get("INDUCTIVA_TASK_RUNNER_IMAGE",
                                   "inductiva/task-runner:main")

FILE_TRACKER_IMAGE = os.environ.get("INDUCTIVA_FILE_TRACKER_IMAGE",
                                    "inductiva/file-tracker:main")

APPTAINER_CONVERTER_IMAGE = os.environ.get(
    "INDUCTIVA_APPTAINER_CONVERTER_IMAGE",
    "inductiva/kutu:apptainer-converter_v0.1.0_dev")

TASK_KILL_MAX_API_REQUESTS = 5

TASK_KILL_RETRY_SLEEP_SEC = 1.0

MAX_CONFIRMATION_LINES = 5

LOADER_COMMAND_PREFIX = "cmd_"

LOADER_IGNORE_PREFIX = "__"

LOADER_HIDE_PREFIX = "_"

INDUCTIVA_GIT_EXAMPLES_URL = ("https://raw.githubusercontent.com/inductiva"
                              "/inductiva/refs/heads/main/inductiva/tests/"
                              "test_simulators/")

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
TMP_DIR = HOME_DIR / "tmp"
TMP_DIR.mkdir(parents=True, exist_ok=True)

TASK_OUTPUT_ZIP = "output.zip"
TASK_INPUT_ZIP = "input.zip"

LOGIN_MESSAGE = "Please login with `inductiva auth login`."
