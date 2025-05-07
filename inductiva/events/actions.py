""" Actions for events"""

from inductiva.client.model.action_email_create import ActionEmailCreate
from inductiva.client.model.action_webhook_create import ActionWebhookCreate


class Action:
    """Base class for actions."""

    def __init__(self) -> None:
        pass

    def get_action(self) -> None:
        raise NotImplementedError("Subclasses must implement this method.")


class EmailNotification(Action):
    """
    Action that sends an email notification.

    Attributes:
        email_address (str): The email address to send the notification to.
    """

    def __init__(self, email_address: str) -> None:
        super().__init__()
        self.email_address = email_address

    def get_action(self) -> ActionEmailCreate:
        """
        Constructs and returns the email action object.

        Returns:
            ActionEmailCreate: Configured email action instance.
        """
        return ActionEmailCreate(
            email_address=self.email_address,
            action_type=ActionEmailCreate.MetaOapg.properties.action_type.EMAIL)


class WebhookNotification(Action):
    """
    Action that sends a webhook notification.

    Attributes:
        webhook_url (str): The webhook url to send the notification to.
    """

    def __init__(self, webhook_url: str) -> None:
        super().__init__()
        self.webhook_url = webhook_url

    def get_action(self) -> ActionEmailCreate:
        """
        Constructs and returns the webhook object.

        Returns:
            ActionWebhookCreate: Configured webhook instance.
        """
        return ActionWebhookCreate(webhook_url=self.webhook_url,
                                   action_type=ActionWebhookCreate.MetaOapg.
                                   properties.action_type.WEBHOOK)
