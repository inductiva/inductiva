"""Module for consuming the streams of a task through a websocket."""
from typing import IO
import logging
import json
import sys

import websocket

import inductiva
from inductiva import constants

logger = logging.getLogger("websocket")
logger.setLevel(logging.INFO)


class TaskStreamConsumer:
    """Consumer of the streams of a running task through a websocket.
    
    Attributes:
        RECONNECT_DELAY_SEC (float): Delay interval when reconnecting in
            seconds.
        PING_INTERVAL_SEC (float): Time interval, in seconds, to automatically
            send a “ping” command  to the websocket (in seconds).
            If set to 0, no ping is sent periodically.
        PING_TIMEOUT_SEC (float): Timeout (in seconds) to close the stream 
            if the pong message is not received.
    """
    RECONNECT_DELAY_SEC = 5.
    PING_INTERVAL_SEC = 15.
    PING_TIMEOUT_SEC = 5.

    def __init__(self,
                 task_id: str,
                 fout: IO = sys.stdout,
                 ferr: IO = sys.stderr):
        """Initialize websocket connection to the task STDOUT & STDERR streams.

        Args:
            task_id (int): ID of the task for which to get the STDOUT & STDERR
                streams.
            fout (IO): I/O streams, such as returned by open(), for the STDOUT.
            ferr (IO): I/O streams, such as returned by open(), for the STDERR.
        """
        query = f'?query={{task_id="{task_id}"}}'
        websocket_url = f"{constants.LOGS_WEBSOCKET_URL}/loki/api/v1/tail"
        self.websocket_url = websocket_url + query
        self.task_id = task_id
        self.fout = fout
        self.ferr = ferr

        # Setup the websocket application to connect to the logs socket.
        self.ws = self.__setup_websocket()

    def __setup_websocket(self):
        """Initialize the websocket app, with the callbacks within the class."""
        return websocket.WebSocketApp(self.websocket_url,
                                      header={"X-API-Key": inductiva.api_key},
                                      on_open=self.__on_open,
                                      on_error=self.__on_error,
                                      on_close=self.__on_close,
                                      on_message=self.__on_message)

    def __on_message(self, unused_ws, message):
        data = json.loads(message)
        streams = data.get("streams")
        for stream in streams:
            print(stream["values"][0][1], file=self.fout)

    def __on_error(self, unused_ws, error):
        if not isinstance(error, KeyboardInterrupt):
            print(error, file=self.ferr)

    def __on_close(self, unused_ws, status_code, message):
        print(f"Closed stream with {status_code}: {message}.", file=self.fout)

    def __on_open(self, unused_ws):
        print(f"Opening socket connection to logs of task {self.task_id} ...",
              file=self.fout)

    def run_forever(self):
        """Stream the STDOUT & STDERR of a task through the websocket."""

        self.ws.run_forever(reconnect=self.RECONNECT_DELAY_SEC,
                            ping_interval=self.PING_INTERVAL_SEC,
                            ping_timeout=self.PING_TIMEOUT_SEC)
