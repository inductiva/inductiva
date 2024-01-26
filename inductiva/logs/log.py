"""Custom logging functions"""
import logging.handlers
import logging
import sys

root_logger = logging.getLogger()


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
    root_logger.error("System encountered the following unhandled exception:\n"
                      "%s\n  Exiting with code 1.", exc_value,
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
        logging.handlers.RotatingFileHandler("rotated.log",
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
