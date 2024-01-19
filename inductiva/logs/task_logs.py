"""Module for streaming the logs of a task through a websocket."""
import sys

import websocket

from inductiva import constants


class TaskLogsStream:
    """Logs the streams of a running task through a websocket."""

    def __init__(self, task_id, fout=sys.stdout, ferr=sys.stderr):
        """Initialize the websocket connection to the task's logs.

        Args:
            task_id (int): ID of the task to print the logs of."""

        self.websocket_url = f"{constants.LOGS_WEBSOCKET}/{task_id}"
        self.task_id = task_id
        self.fout = fout
        self.ferr = ferr
        self.ws = self.setup_websocket()

    def setup_websocket(self):
        """Initialize the websocket app, with the callbacks within the class."""
        return websocket.WebSocketApp(self.websocket_url,
                                      on_open=self.on_open,
                                      on_message=self.on_message,
                                      on_error=self.on_error,
                                      on_close=self.on_close)

    def on_message(self, ws, message):  # pylint: disable=unused-argument
        print(message, file=self.fout)

    def on_error(self, ws, error):  # pylint: disable=unused-argument
        if not isinstance(error, KeyboardInterrupt):
            print(error, file=self.ferr)

    def on_close(self, ws, status_code, message):  # pylint: disable=unused-argument
        print(f"Closed stream with {status_code}: {message}.", file=self.fout)

    def on_open(self, ws):  # pylint: disable=unused-argument
        print("Opening stream of the task logs...", file=self.fout)

    def stream_task_logs(self):
        """Stream the logs of a task through a websocket.
        
        Args:
            task_id (int): ID of the task to print the logs of."""
        self.ws.run_forever(reconnect=5, ping_interval=10, ping_timeout=5)
