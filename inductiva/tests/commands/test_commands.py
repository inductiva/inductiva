import pytest
from pytest import mark

from inductiva import commands


@mark.parametrize("cmd,prompts,expected_json", [
    ("command with prompts", ["y", "y"], {
        "cmd": "command with prompts",
        "prompts": ["y", "y"]
    }),
    ("command without prompts", [], {
        "cmd": "command without prompts",
        "prompts": []
    }),
])
def test_to_json_method(cmd, prompts, expected_json):
    cmd = commands.Command(cmd, *prompts)
    result = cmd.to_json()
    assert result == expected_json


@pytest.mark.parametrize("cmd,prompts,expected_message", [
    ("fails on prompts", ["y", 42], "Prompts must be all strings."),
    (0, ["fails on cmd"], "cmd argument must be a string."),
])
def test_to_json_method_with_invalid_input(cmd, prompts, expected_message):
    with pytest.raises(AssertionError, match=expected_message):
        commands.Command(cmd, *prompts)
