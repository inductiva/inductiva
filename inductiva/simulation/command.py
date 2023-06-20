"""Class for construction of simulator commands."""


class Command():
    """Class for construction of simulator commands."""

    def __init__(self, name: str, **kwargs):
        """Command constructor.
        Args:
            name: The name of the specific method to be executed in 
            the simulation API.
            **kwargs: Additional keyword arguments to be passed to
            the simulation API method."""
        self.name = name
        self.flags = kwargs

    def get_args(self):
        """Return the flags used in the command as a dictionary."""
        all_args = {"name": self.name, **self.flags}

        return all_args
