"""Custom logging functions"""
import os

import logging.handlers
import logging
import sys

import pathlib
import platform

root_logger = logging.getLogger()


def get_logs_file_path():
    system = platform.system()
    logs_name = "inductiva.log"
    if system.lower() == "windows":
        logs_file_path = pathlib.Path.home(
        ) / "AppData" / "inductiva" / logs_name
    elif system.lower() in ["linux", "darwin"]:
        logs_file_path = pathlib.Path.home() / ".inductiva" / logs_name
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
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    root_logger.error(
        "System encountered the following unhandled exception:\n"
        "%s\n  Exiting with code 1.",
        exc_value,
        exc_info=(exc_type, exc_value, exc_traceback))


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
    if not os.path.exists(logs_file_dir):
        os.mkdir(logs_file_dir)

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
    sys.excepthook = handle_uncaught_exception
