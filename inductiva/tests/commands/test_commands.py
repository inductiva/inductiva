"""Tests Command abstraction"""
import pytest
from pytest import mark

from inductiva import commands
from inductiva.commands import Command
from inductiva.commands.mpiconfig import MPIConfig


def test_to_dict__with_prompts():
    cmd = commands.Command("command with prompts", "y", "y")
    assert cmd.to_dict() == {
        "cmd": "command with prompts",
        "prompts": ["y", "y"],
        "mpi_config": None,
        "env": {},
    }


def test_to_dict__without_prompts():
    cmd = commands.Command("command without prompts")
    assert cmd.to_dict() == {
        "cmd": "command without prompts",
        "prompts": [],
        "mpi_config": None,
        "env": {},
    }


@mark.parametrize(
    "command,should_fail",
    [
        ("ls", False),
        ("ls -l", False),
        ("ls | wc", True),
        ("ls > file.txt", True),
        ("sort < file.txt", True),
        ("ls & ls", True),
        ("ls ; ls", True),
        #Some simulators need the char * to be used
        #example "convert -delay 20 -loop 0 potential_*.png scattering.gif"
        ("ls *", False),
        ("ls ?", True),
        ("ls ~", True),
        ("ls $HOME", True),
        ("ls $PATH", True),
    ])
def test_has_special_chars(command, should_fail):
    if should_fail:
        with pytest.raises(ValueError):
            commands.Command(command)
    else:
        # No exception should be raised
        commands.Command(command)


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
        },
        "env": {}
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
        "mpi_config": None,
        "env": {}
    }


def test_to_dict_and_from_dict_roundtrip():
    cfg = MPIConfig("4.1.6", np=4, use_hwthread_cpus=True)

    d = cfg.to_dict()
    assert d["version"] == "4.1.6"
    assert d["options"] == {"np": 4, "use-hwthread-cpus": True}

    # Round-trip check
    cfg2 = MPIConfig.from_dict(d)
    assert cfg2.version == "4.1.6"
    assert cfg2.options == {"np": 4, "use_hwthread_cpus": True}

def test_from_dict_handles_missing_options():
    d = {"version": "1.10.7"}
    cfg = MPIConfig.from_dict(d)
    assert cfg.version == "1.10.7"
    assert cfg.options == {}

def test_dicts_to_commands_basic():
    dicts = [
        {"cmd": "gmx pdb2gmx -f protein.pdb", "prompts": ["y", "y"]},
        {"cmd": "terminate", "prompts": []},
    ]

    cmds = Command.dicts_to_commands(dicts)

    assert len(cmds) == 2
    assert isinstance(cmds[0], Command)
    assert cmds[0].cmd == "gmx pdb2gmx -f protein.pdb"
    assert cmds[0].prompts == ["y", "y"]
    assert cmds[1].cmd == "terminate"
    assert cmds[1].prompts == []

def test_dicts_to_commands_with_mpi_config():
    d = [
        {
            "cmd": "mpirun myprog",
            "prompts": [],
            "mpi_config": {
                "version": "4.1.6",
                "options": {"np": 2, "use-hwthread-cpus": True},
            },
        }
    ]

    cmds = Command.dicts_to_commands(d)
    cmd = cmds[0]

    assert cmd.cmd == "mpirun myprog"
    assert isinstance(cmd.mpi_config, MPIConfig)
    assert cmd.mpi_config.version == "4.1.6"
    assert cmd.mpi_config.options == {"np": 2, "use_hwthread_cpus": True}
