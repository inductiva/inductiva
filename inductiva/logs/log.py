"""Custom logging functions"""
import json
import os

import logging.handlers
import logging
import sys

import pathlib
import platform

from inductiva import constants
from inductiva.client import exceptions
from inductiva.utils import format_utils

root_logger = logging.getLogger()


def get_logs_file_path():
    system = platform.system()
    logs_name = "inductiva.log"
    if system.lower() == "windows":
        logs_file_path = pathlib.Path.home(
        ) / "AppData" / "inductiva" / logs_name
    elif system.lower() in ["linux", "darwin"]:
        logs_file_path = constants.HOME_DIR / logs_name
    else:
        raise RuntimeError(f"Current operating system {system} not supported.")
    return logs_file_path


class NoExceptionFormatter(logging.Formatter):

    def format(self, record):
        record.exc_text = record.exc_info = None
        return super().format(record)


def log_and_exit(logger, level, msg, *args, **kwargs):
    logger.log(level, msg, *args, **kwargs)
    sys.exit(1)


def handle_uncaught_exception(exc_type, exc_value, exc_traceback):
    if _handle_api_exception(exc_type, exc_value, exc_traceback):
        return
    sys.__excepthook__(exc_type, exc_value, exc_traceback)

def _handle_api_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, exceptions.ApiException) and \
        400 <= exc_value.status  < 500:
        detail = json.loads(exc_value.body)["detail"]
        root_logger.error("Error: %s",
                          detail,
                          exc_info=(exc_type, exc_value, exc_traceback))
        return True
    return False

def ipy_handle_uncaught_exception(self,
                                  exc_type,
                                  exc_value,
                                  exc_traceback,
                                  tb_offset=None):
    if _handle_api_exception(exc_type, exc_value, exc_traceback):
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
