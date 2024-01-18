"""Wrapper class for simulator commands."""


class Command(dict):
    """Abstraction class for commands."""

    def __init__(self, cmd, *prompts):
        super().__init__(cmd=cmd, prompts=list(prompts))
