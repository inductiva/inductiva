"""Custom exceptions for tasks"""
import json

from inductiva.client.exceptions import ApiException


class TaskSubmissionException(Exception):
    """
    Custom exception class for handling errors for task submission.

    Attributes:
        exception: The ApiException instance that triggered this exception.
        reason: The reason for the exception, extracted from ApiException.body.
    """

    def __init__(self, exception: ApiException):
        """
        Initializes the TaskSubmissionException.

        Args:
            exception: The original ApiException instance.
        """
        super().__init__()
        self.exception = exception
        self.reason = exception.body

    def __str__(self):
        """
        Returns a user-friendly representation of the error.

        If ApiException.body is a JSON string, extract a detailed error message.
        If parsing the JSON fails or the expected detail key is missing, falls
        back to displaying the ApiException as a string.
        """
        try:
            body_json = json.loads(self.reason)
            message = body_json["detail"]
        except (json.JSONDecodeError, KeyError):
            message = str(self.exception)
        return message
