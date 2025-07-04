""" Triggers for events."""

from pydantic import BaseModel
from inductiva.client.models.trigger_task_type import TriggerTaskType


class Trigger(BaseModel):
    """Base class for triggers."""

    def get_trigger(self) -> None:
        raise NotImplementedError("Subclasses must implement this method.")


class TaskOutputUploaded(Trigger):
    """
    Trigger that is activated when a task output is uploaded.

    Attributes:
        task_id (int): ID of the task to monitor.
    """
    task_id: str

    def get_trigger(self):
        return {
            "task_id": self.task_id,
            "trigger": TriggerTaskType.TASK_OUTPUT_UPLOADED,
            "trigger_type": "task"
        }


class ObserverFileExists(Trigger):
    """
    Trigger that is activated when a file exists.
    
    Attributes:
        task_id (str): ID of the task to monitor.
        file_path (str): Path of the file to monitor."""
    task_id: str
    file_path: str

    def get_trigger(self):
        return {
            "trigger_type": "observer",
            "task_id": self.task_id,
            "trigger": "file_exists_observer",
            "file_path": self.file_path,
        }


class ObserverFileRegex(Trigger):
    """
    Trigger that is activated when a file matches a regex pattern.
    
    Attributes:
        task_id (str): ID of the task to monitor.
        file_path (str): Path of the file to monitor.
        regex (str): Regex pattern to match."""
    task_id: str
    file_path: str
    regex: str

    def get_trigger(self):
        return {
            "trigger_type": "observer",
            "task_id": self.task_id,
            "trigger": "file_regex_observer",
            "file_path": self.file_path,
            "regex": self.regex,
        }
