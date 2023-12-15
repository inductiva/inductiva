"""Mixin for managing commands."""
import io
import json
import os
import typing

from inductiva.utils import templates

Command = typing.NamedTuple("Command", [("cmd", str), ("prompts", list)])


class CommandManager:
    """Class for command management."""

    __commands = None  #pylint: disable=invalid-name

    def add_command(self, command_template, prompts, **render_args):
        """Add commands to the command list."""
        if self.__commands is None:
            self.__commands = []

        command = templates.render_from_string(command_template, **render_args)
        self.__commands.append(Command(command, prompts))

    def add_commands_from_file(self, commands_file, **render_args):
        """Add commands from a file to the command list."""
        inmemory_file = None

        if self.__commands is None:
            self.__commands = []

        if os.path.isfile(commands_file):
            inmemory_file = io.StringIO()
            templates.render_file(source_file=commands_file,
                                  target_file=inmemory_file,
                                  **render_args)

        inmemory_file.seek(0)
        commands = json.load(inmemory_file)
        for command in commands:
            if "prompts" not in command:
                command["prompts"] = []

            self.__commands.append(Command(command["cmd"], Command["prompts"]))
        inmemory_file.close()

    def get_commands(self):
        return self.__commands
