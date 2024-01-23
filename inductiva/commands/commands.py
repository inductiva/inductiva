"""Wrapper class for simulator commands."""


class Command:
    """Abstraction class for commands."""

    def __init__(self, cmd: str, *prompts: str):
        assert isinstance(cmd, str), "cmd argument must be a string."
        for prompt in prompts:
            assert isinstance(prompt, str), "Prompts must be all strings."

        self.cmd = cmd
        self.prompts = list(prompts)

    def to_json(self):
        return self.__dict__
