"""Tests for the templates utils."""
import os
import json

import pytest

from inductiva.utils import templates


@pytest.mark.parametrize("command_template, render_args, command_output",
                         [("grep {{ name }} {{ filename }}", {
                             "name": "time",
                             "filename": "out.log"
                         }, "grep time out.log"),
                          ("echo inductiva is awesome > awesome.txt", {},
                           "echo inductiva is awesome > awesome.txt")])
def test_render_command(command_template, render_args, command_output):
    """Test rendering commands from strings."""

    command = templates.render_from_string(command_template, **render_args)

    assert command == command_output


@pytest.mark.parametrize(
    "file, render_args, output_file, expected_output",
    [(os.path.join(os.path.dirname(__file__), "..", "assets",
                   "commands_template.json.jinja"), {
                       "name": "time"
                   },
      os.path.join(os.path.dirname(__file__), "..", "assets",
                   "command.json"), "echo time logs.txt")])
def test_render_file(file, render_args, output_file, expected_output):
    """Test rendering commands from files."""

    templates.render_file(source_file=file,
                          target_file=output_file,
                          **render_args)

    assert os.path.isfile(output_file)

    with open(output_file, "r", encoding="utf-8") as commands_file:
        commands = json.load(commands_file)
    command = commands[0]
    assert expected_output in command["cmd"]


@pytest.mark.parametrize(
    "template_dir, render_args, target_dir, expected_output",
    [(os.path.join(os.path.dirname(__file__), "..", "assets", "template_dir"), {
        "simulator": "openfoam",
        "simulator_version": "esi"
    }, os.path.join(os.path.dirname(__file__), "..", "assets", "target_dir"), {
        "simulator": "openfoam",
        "simulator_version": "esi"
    })])
def test_render_directory(template_dir, render_args, target_dir,
                          expected_output):
    """Test rendering directories."""

    templates.render_directory(template_dir, target_dir, **render_args)
    target_file = os.path.join(target_dir, "template.json")

    assert os.path.isfile(target_file)

    with open(target_file, "r", encoding="utf-8") as file:
        template_dict = json.load(file)
    assert template_dict == expected_output
