"""Custom logging functions"""
import traceback
import json
import os

import logging.handlers
import logging
import sys

import inductiva
from inductiva import constants
from inductiva.client import exceptions
from inductiva.utils import format_utils

root_logger = logging.getLogger()


def is_cli():
    """Determines if the caller is the CLI"""
    caller = sys.argv[0]
    return caller.endswith(("inductiva", "inductiva.exe"))


def get_logs_file_path():
    return inductiva.constants.LOGS_FILE_PATH


class NoExceptionFormatter(logging.Formatter):

    def format(self, record):
        record.exc_text = record.exc_info = None
        return super().format(record)


def handle_uncaught_exception(exc_type, exc_value, exc_traceback):
    if _handle_api_exception(exc_type,
                             exc_value,
                             exc_traceback,
                             is_notebook=False):
        return
    sys.__excepthook__(exc_type, exc_value, exc_traceback)


def _get_traceback_first_and_last_lines(exc_traceback, n_lines: int,
                                        is_notebook: bool):
    """
    This method gets the first and last N lines of the traceback.
    Args:
        exc_traceback: traceback object.
        n_lines: number of lines to get.
    """
    tb_list = traceback.extract_tb(exc_traceback)

    skip = 1 if is_notebook else 0

    tb_first_list = tb_list[skip:n_lines + skip]
    tb_last_list = tb_list[-n_lines:]

    first_formatted_tb = [
        f"  - {frame.filename} line {frame.lineno}, in {frame.name}"
        for frame in tb_first_list
    ]
    last_formatted_tb = [
        f"  - {frame.filename} line {frame.lineno}, in {frame.name}"
        for frame in tb_last_list
    ]

    result_string = ("\n".join(first_formatted_tb) + "\n" + \
                     "    ...\n" +
                     "\n".join(last_formatted_tb) + "\n")

    return result_string


def _handle_api_exception(exc_type, exc_value, exc_traceback,
                          is_notebook: bool):
    if issubclass(exc_type, exceptions.ApiException) and \
        400 <= exc_value.status  < 500:
        detail = json.loads(exc_value.body)["detail"]

        # Gets the last N lines of the traceback
        formatted_tb = _get_traceback_first_and_last_lines(
            exc_traceback,
            constants.EXCEPTIONS_MAX_TRACEBACK_DEPTH,
            is_notebook=is_notebook)

        if not is_cli():
            detail = (f"{detail}\n  in:\n{formatted_tb}\n"
                      "For more information on this error, "
                      f"check the logs at {get_logs_file_path()}")

        root_logger.error("ERROR: %s",
                          detail,
                          exc_info=(exc_type, exc_value, exc_traceback))
        return True
    if issubclass(exc_type, exceptions.ApiValueError):
        root_logger.error("Error: %s",
                          exc_value,
                          exc_info=(exc_type, exc_value, exc_traceback))
        return True
    if issubclass(exc_type, inductiva.VersionError):
        root_logger.error(exc_value)
        return True

    return False


def ipy_handle_uncaught_exception(self,
                                  exc_type,
                                  exc_value,
                                  exc_traceback,
                                  tb_offset=None):
    if _handle_api_exception(exc_type,
                             exc_value,
                             exc_traceback,
                             is_notebook=True):
        return
    self.showtraceback((exc_type, exc_value, exc_traceback),
                       tb_offset=tb_offset)


def setup(level=logging.INFO):
    formatters = [
        logging.Formatter(
            fmt=
            "%(asctime)s|%(levelname)s|%(name)s|%(filename)s:%(lineno)"\
            "d %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"),
        NoExceptionFormatter(fmt="%(message)s")
    ]

    logs_file_path = get_logs_file_path()
    logs_file_dir = logs_file_path.parent
    os.makedirs(logs_file_dir, exist_ok=True)

    handlers = [
        logging.handlers.RotatingFileHandler(logs_file_path,
                                             encoding="utf8",
                                             maxBytes=1e6,
                                             backupCount=10),
        logging.StreamHandler()
    ]

    for handler, formatter in zip(handlers, formatters):
        handler.setFormatter(formatter)
        root_logger.addHandler(handler)

    root_logger.setLevel(level)

    api_traceback = format_utils.getenv_bool("INDUCTIVA_DEBUG_API_TRACEBACK",
                                             False)
    if api_traceback:
        return

    has_ipython = False
    try:
        ipython = __import__("IPython")
        has_ipython = True
    except ImportError:
        pass

    if has_ipython and (ip := ipython.get_ipython()):
        ip.set_custom_exc((Exception,), ipy_handle_uncaught_exception)
    else:
        sys.excepthook = handle_uncaught_exception
