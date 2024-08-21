"""Wrapper class for simulator commands."""
from inductiva import types
from .mpiconfig import MPIConfig


class Command:
    """Abstraction class for commands."""

    def __init__(self, cmd: str, *prompts: str, mpi_config: MPIConfig = None):
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

        if mpi_config and not isinstance(mpi_config, MPIConfig):
            raise TypeError("'mpi_config' must be an instance of MPIConfig.")

        self.cmd = cmd
        self.prompts = list(prompts)
        self.mpi_config = mpi_config

    def to_dict(self):
        """
        Convert this command to a dictionary representation.
        
        Returns:
        A dictionary with {"cmd": self.cmd, "prompts": self.prompts}

        Example:
        >>> cmd = Command("gmx pdb2gmx -f protein.pdb", "y", "y")
        >>> cmd.to_dict()
        >>> {"cmd": "gmx pdb2gmx -f protein.pdb", prompts: ["y", "y"]}
        """
        ret_dict = self.__dict__.copy()

        if ret_dict["mpi_config"] is not None:
            ret_dict["mpi_config"] = ret_dict["mpi_config"].to_dict()

        return ret_dict

    @staticmethod
    def commands_to_dicts(commands: types.Commands):
        """
        Convert the given commands to a list of dictionaries.

        Example:
        >>> commands = [Command("run", "y"), "terminate"]
        >>> simulator.commands_to_dict(commands)
        >>> [{"cmd": "run", "prompts": ["y"]},
        ...  {"cmd": "terminate", "prompts": []}]
        """
        return [
            cmd.to_dict()
            if isinstance(cmd, Command) else Command(cmd).to_dict()
            for cmd in commands
        ]
