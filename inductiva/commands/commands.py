"""Wrapper class for simulator commands."""


class Command:
    """Abstraction class for commands."""

    def __init__(self, cmd: str, *prompts: str):
        if not isinstance(cmd, str):
            raise ValueError("cmd argument must be a string.")

        if not all(isinstance(prompt, str) for prompt in prompts):
            raise ValueError("Prompts must be all strings.")

        self.cmd = cmd
        self.prompts = list(prompts)

    def to_dict(self):
        return self.__dict__.copy()
