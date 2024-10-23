"""Module for consuming the streams of a task through a websocket.

Task logs are output in the following format:

$ induction logs <task_id>
<log line 1>                                        log content
<log line 2>
  ...
<log line N>
                                                    `SPACING` black lines
 ● 0:00:02 Streaming output from task <task_id>     connection status footer
"""
from threading import Thread, RLock
from urllib.parse import urlencode
from datetime import timedelta
from typing import IO
import logging
import sched
import json
import sys
import time
import os

import websocket

import inductiva
from inductiva import constants

logger = logging.getLogger("websocket")
logger.setLevel(logging.CRITICAL)

SPACING = 1
ESC = "\u001B["
RESET = f"{ESC}0m"
EMPHASIS = f"{ESC}2m"
FAILURE = f"{ESC}31m"
SUCCESS = f"{ESC}32m"
PENDING = f"{ESC}33m"
CLEAR_LINE = f"{ESC}0K"
UP = f"{ESC}{SPACING}A"
DOWN = f"{ESC}{SPACING}B"

ANSI_ENABLED = os.getenv("ANSI_ENABLED", "1") == "1"


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
    TICK_INTERVAL_SEC = 1.
    END_OF_STREAM = "<<end_of_stream>>"
    ALLOWED_IO_STREAMS = {
        "std_err": "standard error",
        "std_out": "standard output"
    }

    def __init__(self,
                 task_id: str,
                 fout: IO = sys.stdout,
                 ferr: IO = sys.stderr,
                 io_stream: str = None,
                 no_color: bool = False):
        """Initialize websocket connection to the task STDOUT & STDERR streams.

        Args:
            task_id (int): ID of the task for which to get the STDOUT & STDERR
                streams.
            fout (IO): I/O streams, such as returned by open(), for the STDOUT.
            ferr (IO): I/O streams, such as returned by open(), for the STDERR.
            io_stream (str): The type of stream to consume. It can be either
                "std_out" or "std_err". If not provided, both are consumed.
            no_color (bool): Disable ANSI color codes in the output. Default is
                False.
        """
        self.task_id = task_id
        self.fout = fout
        self.ferr = ferr

        # Check if the file handles are interactive and if ANSI is enabled.
        # Note that there is no point in using ANSI if ferr is not interactive
        self._ferr_isatty = ANSI_ENABLED and self.ferr.isatty()
        self._fout_isatty = ANSI_ENABLED and self.fout.isatty()
        self._fout_isatty &= self._ferr_isatty
        self._message_formatter = self._get_message_formatter()
        self._status_formatter = self._get_status_formatter()

        # housekeeping variables
        self._conn_start_time = None
        self._conn_alive_sec = 0.0
        self._conn_opened = False
        self._prev_status = None
        self._conn_count = 0
        self._lock = RLock()

        self._last_message_timestamp = None
        self._keyboard_interrupt = False
        self._end_of_stream = False

        if io_stream and io_stream not in self.ALLOWED_IO_STREAMS:
            raise ValueError(
                f"Invalid io_stream: {io_stream}\n"
                f"Allowed values are: {list(self.ALLOWED_IO_STREAMS.keys())}.")
        self.io_stream = io_stream
        self.no_color = no_color

    def __on_message(self, ws: websocket.WebSocketApp, message: str):
        data = json.loads(message)
        streams = data.get("streams")
        for stream in streams:
            self._last_message_timestamp, msg = stream["values"][0]

            if msg == self.END_OF_STREAM:
                self._end_of_stream = True
                ws.close()
                return
            # Only writes in red if we are receiving stdout and stderr
            # and the output is stderr
            is_stderr = (self.io_stream is None and
                         stream["stream"]["io_type"] == "IOTypes.STD_ERR")
            self._write_message(msg.rstrip(),
                                is_stderr=is_stderr and not self.no_color)

    def _get_message_formatter(self):
        """Get the message formatter based on the output file
        handle interactivity.
        """
        if self._fout_isatty:
            return self._ansi_message_formatter
        return self._plain_message_formatter

    def _get_status_formatter(self):
        """Get the message formatter based on the error file
        handle interactivity.
        """
        if self._ferr_isatty:
            return self._ansi_status_formatter
        return self._plain_status_formatter

    def _ansi_message_formatter(self, msg, is_stderr=False):
        """Fancy message formatter for formatted message output using
        ANSI escape codes.
        """
        # 1. move up SPACING lines
        # 2. move to beginning of line and clear line
        # 3. write message
        # 4. add new line
        # 5. move down SPACING lines
        # 6. move to beginning of line and clear line
        # 7. add new blank line
        if is_stderr:
            m = f"{UP}\r{CLEAR_LINE}{FAILURE}{msg}{RESET}{DOWN}\r{CLEAR_LINE}\n"
        else:
            m = f"{UP}\r{CLEAR_LINE}{msg}{DOWN}\r{CLEAR_LINE}\n"
        return m

    def _plain_message_formatter(self, msg, **_):
        """Simple message formatter for plain message output."""
        return msg + "\n"

    def _ansi_status_formatter(self, msg, status):
        """Fancy message formatter for formatted status message output using
        ANSI escape codes.
        """
        if self._conn_start_time is None:
            elapsed = ""
        else:
            alive = timedelta(seconds=int(self._conn_alive_sec))
            elapsed = f"{str(alive)} "

        # 1. go to beginning of line
        # 2. clear line
        # 3. set color for bullet
        #   3.1 write bullet icon
        #   3.2 reset color
        # 4. set emphasis formatter
        #   4.1. write timer message (if any)
        #   4.2. write status message
        # 5. reset everything
        s = f"\r{CLEAR_LINE}{status} ● {RESET}{EMPHASIS}{elapsed}{msg}{RESET}"
        return s

    # _ to make it have the same number or arguments as the other formatters
    # In this case, the second argument is not used.
    def _plain_status_formatter(self, msg, _):
        """Simple message formatter for plain status message output."""
        return msg + "\n"

    def _write_message(self, msg, is_stderr=False):
        """Write the message to the output file handle

        Args:
            msg (str): message to be written to output.
        """
        s = self._message_formatter(msg, is_stderr=is_stderr)
        self.fout.write(s)
        self.fout.flush()
        self._redraw_status()

    def _write_status(self, msg, status):
        """Write the status message to the error file handle

        Args:
            msg (str): Status message.
            status (str): Status formatting code.
        """
        s = self._status_formatter(msg, status)
        self.ferr.write(s)
        self.ferr.flush()

    def _update_status(self, msg, status):
        """Update the connection status message.

        Args:
            msg (str): Status message.
            status (str): Status formatting code.
        """
        with self._lock:
            self._prev_status = (msg, status)
            self._write_status(msg, status)

    def _redraw_status(self):
        """Redraw the connection status message."""
        if self._prev_status and self._ferr_isatty:
            self._update_status(*self._prev_status)

    def __on_error(self, ws: websocket.WebSocketApp, error: Exception):
        """Callback for websocket error event."""
        if isinstance(error, KeyboardInterrupt):
            self._keyboard_interrupt = True
            ws.close()
        else:
            self._update_status(str(error), FAILURE)
            time.sleep(2)

    def __on_close(self, unused_ws, close_status_code, close_message):
        """Callback for websocket close event."""

        if close_status_code or close_message:
            close_msg = f" (code: {close_status_code}, msg: {close_message})"
        else:
            close_msg = ""

        if self._conn_opened:
            if self._end_of_stream:
                msg = f"Received end of stream from task {self.task_id}"
                status = SUCCESS
            else:
                msg = f"Disconnected from task {self.task_id}"
                status = SUCCESS if self._keyboard_interrupt else FAILURE
            msg = f"{msg}{close_msg}."
        else:
            msg = f"Failed to connect to task {self.task_id}{close_msg}."
            status = FAILURE

        self._update_status(msg, status)

    def __on_open(self, unused_ws):
        """Callback for websocket open event."""
        self._conn_start_time = time.time()
        self._conn_opened = True
        self._conn_count += 1

        output = self.ALLOWED_IO_STREAMS.get(self.io_stream, "output")
        msg = f"Streaming {output} from task {self.task_id}"
        if ANSI_ENABLED and not self._fout_isatty:
            msg += " (redirected)"

        self._update_status(msg, SUCCESS)

    def __init_status(self):
        """Allocate sufficient lines for the status message."""
        if self._ferr_isatty:
            self.ferr.write(SPACING * "\n")
            self.ferr.flush()

    def _run_websocket(self):
        """Run the websocket application.
        If connection is lost, reconnect after a delay
        and update the start query to get the logs from the last
        received timestamp."""

        endpoint = constants.LOGS_WEBSOCKET_URL + "/loki/api/v1/tail?"
        if self.io_stream:
            io_type = f'io_type="IOTypes.{self.io_stream.upper()}"'
            params = {
                "query": f'{{{io_type},task_id="{self.task_id}"}}',
                "limit": 500
            }
        else:
            params = {"query": f'{{task_id="{self.task_id}"}}', "limit": 500}
        while True:
            url = endpoint + urlencode(params)
            ws = websocket.WebSocketApp(
                url,
                header={"X-API-Key": inductiva.get_api_key()},
                on_open=self.__on_open,
                on_error=self.__on_error,
                on_close=self.__on_close,
                on_message=self.__on_message)

            ws.run_forever(reconnect=0,
                           ping_interval=self.PING_INTERVAL_SEC,
                           ping_timeout=self.PING_TIMEOUT_SEC)

            if self._end_of_stream or self._keyboard_interrupt:
                break

            if self._last_message_timestamp is not None:
                params["start"] = int(self._last_message_timestamp) + 1
                params["limit"] = 10000

            time.sleep(self.RECONNECT_DELAY_SEC)

    def run_forever(self):
        """Stream the STDOUT & STDERR of a task through the websocket."""

        msg = f"Opening connection to task {self.task_id}..."
        self.__init_status()
        self._update_status(msg, PENDING)
        time.sleep(2)

        if scheduler := self._get_ontick_scheduler():
            thread = Thread(target=scheduler.run)
            thread.daemon = True
            thread.start()

        self._run_websocket()

    def _get_ontick_scheduler(self):
        """Get a scheduler for the connection alive ticker.

        Gets a scheduler that runs every second to update the time elapsed
        since the connection was opened.
        """
        if not self._ferr_isatty:
            return None

        scheduler = sched.scheduler(time.time, time.sleep)
        scheduler.enter(self.TICK_INTERVAL_SEC, 1, self._ontick, (scheduler,))
        return scheduler

    def _ontick(self, scheduler):
        """Update the connection alive time every second.

        Updates the time elapsed since the connection was opened, redraws
        the status message and schedules the next update.
        """
        if self._conn_start_time is not None:
            self._conn_alive_sec = time.time() - self._conn_start_time
            self._redraw_status()

        # schedule the next update
        scheduler.enter(self.TICK_INTERVAL_SEC, 1, self._ontick, (scheduler,))
