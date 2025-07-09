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


class TriggerCreditsInfo:
    """Model for the credits trigger."""

    def __init__(self, trigger_type: str, trigger: str):
        self.trigger_type = trigger_type
        self.trigger = trigger

    def __str__(self):
        return (f"TriggerCreditsInfo(type={self.trigger_type}, "
                f"trigger={self.trigger})")

    __repr__ = __str__


class TriggerFileExistsObserverInfo:
    """Model for the observer trigger."""

    def __init__(self, trigger_type: str, observer_id: uuid.UUID,
                 observer_type: str, task_id: str, file_path: str):
        self.trigger_type = trigger_type
        self.observer_id = observer_id
        self.observer_type = observer_type
        self.task_id = task_id
        self.file_path = file_path

    def __str__(self):
        return (f"TriggerFileExistsObserverInfo("
                f"task_id={self.task_id}, "
                f"file_path={self.file_path})")

    __repr__ = __str__


class TriggerRegexObserverInfo:
    """Model for the observer trigger."""

    def __init__(self, trigger_type: str, observer_id: uuid.UUID,
                 observer_type: str, task_id: str, file_path: str, regex: str):
        self.trigger_type = trigger_type
        self.observer_id = observer_id
        self.observer_type = observer_type
        self.task_id = task_id
        self.file_path = file_path
        self.regex = regex

    def __str__(self):
        return (f"TriggerRegexObserverInfo("
                f"task_id={self.task_id}, "
                f"file_path={self.file_path}), "
                f"regex={self.regex})")

    __repr__ = __str__


ActionInfoUnion = Union[ActionEmailInfo, ActionWebhookInfo]
TriggerInfoUnion = Union[TriggerTaskInfo, TriggerMachineGroupInfo,
                         TriggerCreditsInfo, TriggerFileExistsObserverInfo,
                         TriggerRegexObserverInfo]


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
    elif trigger_type == "credits":
        trigger = TriggerCreditsInfo(**trigger_data)
    elif trigger_type == "observer":
        observer_type = trigger_data.get("observer_type")
        if observer_type == "file_exists_observer":
            trigger = TriggerFileExistsObserverInfo(**trigger_data)
        elif observer_type == "file_regex_observer":
            trigger = TriggerRegexObserverInfo(**trigger_data)
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
