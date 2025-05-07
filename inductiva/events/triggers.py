""" Triggers for events."""

from inductiva.client.model.trigger_task_create import TriggerTaskCreate
from inductiva.client.model.trigger_task_type import TriggerTaskType


class Trigger:
    """Base class for triggers."""

    def __init__(self) -> None:
        pass

    def get_trigger(self) -> None:
        raise NotImplementedError("Subclasses must implement this method.")


class TaskOutputUploaded(Trigger):
    """
    Trigger that is activated when a task output is uploaded.
    
    Attributes:
        task_id (int): ID of the task to monitor.
    """

    def __init__(self, task_id: str) -> None:
        super().__init__()
        self.task_id = task_id

    def get_trigger(self) -> TriggerTaskCreate:
        """
        Constructs and returns the trigger object for output upload.

        Returns:
            TriggerTaskCreate: Configured trigger instance.
        """
        return TriggerTaskCreate(task_id=self.task_id,
                                 trigger=TriggerTaskType.OUTPUT_UPLOADED,
                                 trigger_type=TriggerTaskCreate.MetaOapg.
                                 properties.trigger_type.TASK)
