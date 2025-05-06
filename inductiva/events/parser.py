""" Parser for the events."""
import uuid
from typing import Union, Any


class ActionEmailInfo:
    """Model for the email action."""

    def __init__(self, action_type: str, email_address: str):
        self.action_type = action_type
        self.email_address = email_address

    def __str__(self):
        return (f"ActionEmailInfo(type={self.action_type}, "
                f"email={self.email_address})")

    __repr__ = __str__


class ActionWebhookInfo:
    """Model for the webhook action."""

    def __init__(self, action_type: str, webhook_url: str):
        self.action_type = action_type
        self.webhook_url = webhook_url

    def __str__(self):
        return (f"ActionWebhookInfo(type={self.action_type}, "
                f"webhook={self.webhook_url})")

    __repr__ = __str__


class TriggerTaskInfo:
    """ Model for the task trigger."""

    def __init__(self, trigger_type: str, task_id: str):
        self.trigger_type = trigger_type
        self.task_id = task_id

    def __str__(self):
        return (f"TriggerTaskInfo(type={self.trigger_type}, "
                f"task_id={self.task_id})")

    __repr__ = __str__


class TriggerMachineGroupInfo:
    """Model for the machine group trigger."""

    def __init__(self, trigger_type: str, machine_group_id: int):
        self.trigger_type = trigger_type
        self.machine_group_id = machine_group_id

    def __str__(self):
        return (f"TriggerMachineGroupInfo(type={self.trigger_type}, "
                f"machine_group_id={self.machine_group_id})")

    __repr__ = __str__


ActionInfoUnion = Union[ActionEmailInfo, ActionWebhookInfo]
TriggerInfoUnion = Union[TriggerTaskInfo, TriggerMachineGroupInfo]


class EventInfo:
    """Model for the event."""

    def __init__(self, event_id: str, trigger: TriggerInfoUnion,
                 action: ActionInfoUnion):
        self.event_id = uuid.UUID(event_id)
        self.trigger = trigger
        self.action = action

    def __str__(self):
        return (f"EventInfo(\n"
                f"  event_id={self.event_id},\n"
                f"  trigger={self.trigger},\n"
                f"  action={self.action}\n"
                f")")

    __repr__ = __str__


def parse_event_info(data: dict[str, Any]) -> EventInfo:
    """Parse the raw API response into an EventInfo object."""

    # Parsing the trigger
    trigger_data = data["trigger"]
    trigger_type = trigger_data.get("trigger_type")

    if trigger_type == "task":
        trigger = TriggerTaskInfo(**trigger_data)
    elif trigger_type == "machine_group":
        trigger = TriggerMachineGroupInfo(**trigger_data)
    else:
        raise ValueError(f"Unknown trigger_type: {trigger_type}")

    # Parsing the action
    action_data = data["action"]
    action_type = action_data.get("action_type")

    if action_type == "email":
        action = ActionEmailInfo(**action_data)
    elif action_type == "webhook":
        action = ActionWebhookInfo(**action_data)
    else:
        raise ValueError(f"Unknown action_type: {action_type}")

    # Return the EventInfo model
    return EventInfo(
        event_id=data["event_id"],
        trigger=trigger,
        action=action,
    )
