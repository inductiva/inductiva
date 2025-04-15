from inductiva.client.model.action_email_create import ActionEmailCreate


class Action:
    """Base class for actions."""

    def __init__(self) -> None:
        pass

    def get_action() -> None:
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
