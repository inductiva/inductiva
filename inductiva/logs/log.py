"""Custom logging functions"""
import os

import logging.handlers
import logging
import sys

import pathlib
import platform

root_logger = logging.getLogger()

system = platform.system()
logs_name = "inductiva.log"

if system == "Windows":
    logs_dir = pathlib.Path.home() / "AppData" / "inductiva"
    LOGS_FILE = logs_dir / logs_name
elif system in ["Linux", "Darwin"]:
    logs_dir = pathlib.Path.home() / ".inductiva"
    LOGS_FILE = logs_dir / logs_name
else:
    raise RuntimeError(f"Current operating system {system} not supported.")

if not os.path.exists(logs_dir):
    os.mkdir(logs_dir)


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

    handlers = [
        logging.handlers.RotatingFileHandler(LOGS_FILE,
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
