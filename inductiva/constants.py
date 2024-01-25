"""Constants that can be set for the Inductiva client."""
import os

LOGS_WEBSOCKET_URL = os.environ.get("INDUCTIVA_TASK_LOGS_URL",
                                    "wss://logs.inductiva.ai/loki/api/v1/tail")

DEFAULT_QUEUE_MACHINE_TYPE = "c2-standard-4"
