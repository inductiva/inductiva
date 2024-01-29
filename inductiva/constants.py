"""Constants that can be set for the Inductiva client."""
import os

LOGS_WEBSOCKET_URL = os.environ.get("INDUCTIVA_TASK_LOGS_URL",
                                    "ws://127.0.0.1:5000/logs")

DEFAULT_QUEUE_MACHINE_TYPE = "c2-standard-4"

MAX_API_RETRIES = 5
