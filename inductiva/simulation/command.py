"""Classe for construction of simulator commands."""
from abc import ABC


class Command(ABC):
    """Abstract class for construction of simulator commands."""

    def __init__(self, **kwargs):
        """Command constructor.
        Args:
            method_name: The name of the method to be executed in 
            the simulation API.
            **kwargs: Additional keyword arguments to be passed to
            the simulation API method."""

        self.kwargs = kwargs

    def get_args(self):
        """Return the command as a dictionary."""
        return self.kwargs
