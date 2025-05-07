"""Inductiva Exceptions"""


class InductivaException(Exception):
    """Base class for all Inductiva exceptions."""

    def __init__(self, message, code=None):
        super().__init__(message)
        self.code = code
