"""Custom logging functions"""
from pathlib import Path
import linecache
import traceback
import json
import os

import logging.handlers
import logging
import sys

import inductiva
from inductiva.client import exceptions
from inductiva.utils import format_utils
from inductiva.utils import InductivaException

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


def _format_traceback_to_string(exc_traceback, is_notebook: bool):
    """
    This method formats the traceback into a nice string.
    Args:
        exc_traceback: traceback object.
    """
    tb_list = traceback.extract_tb(exc_traceback)

    skip = 1 if is_notebook else 0

    tb_list = tb_list[skip:]

    formatted_tb = [
        f"{frame.filename} line {frame.lineno}" for frame in tb_list
    ]

    result_string = ("\n".join(formatted_tb))

    return result_string


def _get_traceback_first_line(exc_traceback, is_notebook: bool):
    """Gets the first line of the traceback.

    The first line has the most relevant information, that is, where
    the error happened.

    returns:
        - The first line of the traceback
    """
    formatted_tb = _format_traceback_to_string(exc_traceback,
                                               is_notebook=is_notebook)

    return formatted_tb.splitlines()[0]


def _log_error(detail, exc_type, exc_value, exc_traceback):
    """Logs the error with the appropriate details."""
    root_logger.error("ERROR: %s",
                      detail,
                      exc_info=(exc_type, exc_value, exc_traceback))


def _format_detail_message(detail, context, formatted_tb):
    return (f"{detail}\n\n  In {formatted_tb}\n\n"
            # f"{context}\n\n"
            "For more information on this error, "
            f"check the logs at {get_logs_file_path()}")


def _is_exception_from_inductiva(exc_tb):
    # Get the root path of your package
    inductiva_path = Path(inductiva.__file__).resolve().parent
    # Walk through the traceback to find if any frame belongs to inductiva
    while exc_tb:
        frame_path = Path(exc_tb.tb_frame.f_code.co_filename).resolve()
        # Check if the frame path is a subpath of the inductiva path
        if os.path.commonpath([str(frame_path),
                               str(inductiva_path)]) == str(inductiva_path):
            return True
        exc_tb = exc_tb.tb_next
    return False


def _get_traceback_context(exc_traceback, depth, context):
    current_depth = 0
    result_string = ""
    while exc_traceback is not None:
        if current_depth != depth:
            current_depth += 1
            exc_traceback = exc_traceback.tb_next
            continue

        frame = exc_traceback.tb_frame
        lineno = exc_traceback.tb_lineno
        filename = frame.f_code.co_filename

        start = max(1, lineno - context)
        end = lineno + context + 1

        for i in range(start, end):
            code_line = linecache.getline(filename, i)
            if code_line:
                pointer = " -->" if i == lineno else "    "
                result_string += f"{pointer} {i}: {code_line.rstrip()}\n"

        exc_traceback = exc_traceback.tb_next
        current_depth += 1
    return result_string


def _handle_api_exception(exc_type, exc_value, exc_traceback,
                          is_notebook: bool):

    # In notebooks we ignore the first frame of the traceback
    context_depth = 1 if is_notebook else 0

    if issubclass(exc_type, exceptions.ApiException) and \
        400 <= exc_value.status  < 500:
        detail = json.loads(exc_value.body)["detail"]

        tb_first_line = _get_traceback_first_line(exc_traceback, is_notebook)

        context = _get_traceback_context(exc_traceback, context_depth, 1)
        if not is_cli():
            detail = _format_detail_message(detail, context, tb_first_line)

        _log_error(detail, exc_type, exc_value, exc_traceback)
        return True
    if issubclass(exc_type, exceptions.ApiValueError):
        _log_error(exc_value, exc_type, exc_value, exc_traceback)

        return True
    if issubclass(exc_type, inductiva.VersionError):
        root_logger.error(exc_value)
        return True

    if issubclass(exc_type, InductivaException):
        detail = str(exc_value)

        tb_first_line = _get_traceback_first_line(exc_traceback, is_notebook)
        context = _get_traceback_context(exc_traceback, context_depth, 1)
        if not is_cli():
            detail = _format_detail_message(detail, context, tb_first_line)

        _log_error(detail, exc_type, exc_value, exc_traceback)

        #Inductiva exit would be called here
        #Would make some cleanup if needed
        #Like turning off some machines etc.

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
