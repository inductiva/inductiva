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
