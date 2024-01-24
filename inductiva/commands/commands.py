"""Wrapper class for simulator commands."""


class Command:
    """Abstraction class for commands."""

    def __init__(self, cmd: str, *prompts: str):
        """
        Args:
          cmd: A string with the command, e.g, 'gmx pdb2gmx -f protein.pdb'
          *promts: Strings to be used as answers to prompts, e.g, 'y'

        Example:
        >>> cmd = Command("gmx pdb2gmx -f protein.pdb", "y", "y")
        """
        if not isinstance(cmd, str):
            raise ValueError("cmd argument must be a string.")

        if not all(isinstance(prompt, str) for prompt in prompts):
            raise ValueError("Prompts must be all strings.")

        self.cmd = cmd
        self.prompts = list(prompts)

    def to_dict(self):
        """
        Returns:
        A dictionary with {"cmd": self.cmd, "prompts": self.prompts}

        Example:
        >>> cmd = Command("gmx pdb2gmx -f protein.pdb", "y", "y")
        >>> cmd.to_dict()
        >>> {"cmd": "gmx pdb2gmx -f protein.pdb", prompts: ["y", "y"]}
        """
        return self.__dict__.copy()
