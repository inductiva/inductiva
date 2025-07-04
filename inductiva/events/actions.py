""" Actions for events"""

from pydantic import BaseModel


class Action(BaseModel):
    """Base class for actions."""

    def get_action(self) -> None:
        raise NotImplementedError("Subclasses must implement this method.")


class EmailNotification(Action):
    """
    Action that sends an email notification.

    Attributes:
        email_address (str): The email address to send the notification to.
    """
    email_address: str

    def get_action(self):
        return {"email_address": self.email_address, "action_type": "email"}


class WebhookNotification(Action):
    """
    Action that sends a webhook notification.

    Attributes:
        webhook_url (str): The webhook url to send the notification to.
    """
    webhook_url: str

    def get_action(self):
        return {"webhook_url": self.webhook_url, "action_type": "webhook"}
