"""Constants that can be set for the Inductiva client."""
import os

LOGS_WEBSOCKET = os.environ.get("INDUCTIVA_TASK_LOGS_URL",
                                "ws://127.0.0.1:5000/logs")
