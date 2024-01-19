"""Wrapper class for simulator commands."""


class Command:
    """Abstraction class for commands."""

    def __init__(self, cmd, *prompts):
        self.cmd = cmd
        self.prompts = list(prompts)

    def to_json(self):
        return self.__dict__
