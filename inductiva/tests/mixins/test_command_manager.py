import os

import pytest

from inductiva import mixins


@pytest.mark.parametrize("command_template, render_args, command_output", [
    pytest.param("grep {{ name }} {{ filename }}",
                 dict(name="time", filename="out.log"), "grep time out.log"),
    pytest.param("echo inductiva is awesome > awesome.txt", dict(),
                 "echo inductiva is awesome > awesome.txt")
],
                         ids=["grep", "echo"])
def test_add_command(command_template, render_args, command_output):
    """Test the command manager."""

    cmd_manager = mixins.command_manager.CommandManager()
    cmd_manager.add_command(command_template, prompts=[], **render_args)

    commands = cmd_manager.get_commands()

    # Check the latest command added to the list
    # Parameteric tests seem to be run in one go.
    assert command_output == commands[0].cmd


@pytest.mark.parametrize(
    "commands_file, render_args, command_output",
    [(os.path.join(os.path.dirname(__file__),
                   "..", "assets", "commands_template.json.jinja"),
      dict(name="time"), "echo time logs.txt")])
def test_add_commands_from_file(commands_file, render_args, command_output):

    cmd_manager = mixins.command_manager.CommandManager()
    cmd_manager.add_commands_from_file(commands_file, **render_args)
    commands = cmd_manager.get_commands()

    # Check the latest command added to the list
    # Parameteric tests seem to be run in one go.
    assert command_output == commands[0].cmd
