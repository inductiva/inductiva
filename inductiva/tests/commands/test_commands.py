"""Tests Command abstraction"""
import pytest
from pytest import mark

from inductiva import commands


def test_to_dict__with_prompts():
    cmd = commands.Command("command with prompts", "y", "y")
    assert cmd.to_dict() == {
        "cmd": "command with prompts",
        "prompts": ["y", "y"]
    }


def test_to_dict__without_prompts():
    cmd = commands.Command("command without prompts")
    assert cmd.to_dict() == {"cmd": "command without prompts", "prompts": []}


@pytest.mark.parametrize("cmd,prompts,expected_message", [
    ("fails on prompts", ["y", 42], "Prompts must be all strings."),
    (0, ["fails on cmd"], "cmd argument must be a string."),
])
def test_to_dict_method_with_invalid_input(cmd, prompts, expected_message):
    with pytest.raises(ValueError, match=expected_message):
        commands.Command(cmd, *prompts)


def test_to_dict_does_not_return_object():
    x = commands.Command("command", "y", "n")
    assert id(x.to_dict()) != id(x.__dict__)
