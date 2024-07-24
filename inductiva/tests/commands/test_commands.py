"""Tests Command abstraction"""
import pytest
from pytest import mark

from inductiva import commands
from inductiva.commands.mpiconfig import MPIConfig


def test_to_dict__with_prompts():
    cmd = commands.Command("command with prompts", "y", "y")
    assert cmd.to_dict() == {
        "cmd": "command with prompts",
        "prompts": ["y", "y"],
        "mpi_config": None
    }


def test_to_dict__without_prompts():
    cmd = commands.Command("command without prompts")
    assert cmd.to_dict() == {
        "cmd": "command without prompts",
        "prompts": [],
        "mpi_config": None
    }


@mark.parametrize("cmd,prompts,expected_message", [
    ("fails on prompts", ["y", 42], "Prompts must be all strings."),
    (0, ["fails on cmd"], "cmd argument must be a string."),
])
def test_to_dict_method_with_invalid_input(cmd, prompts, expected_message):
    with pytest.raises(ValueError, match=expected_message):
        commands.Command(cmd, *prompts)


def test_to_dict_does_not_return_object():
    x = commands.Command("command", "y", "n")
    assert id(x.to_dict()) != id(x.__dict__)


def test_to_dict__with_mpi_config():
    config = MPIConfig("1.2.3", np=4)
    cmd = commands.Command("command with prompts", "y", "y", mpi_config=config)
    assert cmd.to_dict() == {
        "cmd": "command with prompts",
        "prompts": ["y", "y"],
        "mpi_config": {
            "version": "1.2.3",
            "options": {
                "np": 4
            }
        }
    }


def test_to_dict__mpi_config_is_not_instance_of_mpiconfig():
    with pytest.raises(TypeError, match="instance of MPIConfig"):
        commands.Command("command with prompts",
                         "y",
                         "y",
                         mpi_config="not an instance of MPIConfig")


def test_to_dict__mpi_config_is_none():
    cmd = commands.Command("command with prompts", "y", "y", mpi_config=None)
    assert cmd.to_dict() == {
        "cmd": "command with prompts",
        "prompts": ["y", "y"],
        "mpi_config": None
    }
