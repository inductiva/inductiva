"""Test for the CLI commands."""
import subprocess
import pytest
import os

#read INDUCTIVA_API_KEY from os env vars
API_KEY = os.environ["INDUCTIVA_API_KEY"]

# ([command], user_input)
CLI_COMMANDS = [
    (["inductiva"], None),
    (["inductiva", "auth"], None),
    (["inductiva", "auth", "login"], API_KEY),
    (["inductiva", "auth", "logout"], None),

    #cant test inductiva logs
    #(["inductiva", "logs"], None, True),
    (["inductiva", "projects"], None),
    (["inductiva", "projects", "ls"], None),
    (["inductiva", "quotas"], None),
    (["inductiva", "quotas", "list"], None),
    (["inductiva", "resources"], None),
    ([
        "inductiva",
        "resources",
        "start",
        "c2-standard-4",
    ], None),
    (["inductiva", "resources", "cost", "c2-standard-4"], None),
    #cant test inductiva resources info
    (["inductiva", "resources", "available"], None),
    (["inductiva", "resources", "ls"], None),
    (["inductiva", "resources", "terminate", "--all"], "yes"),
    (["inductiva", "simulators"], None),
    (["inductiva", "simulators", "list"], None),
    (["inductiva", "storage"], None),
    (["inductiva", "storage", "ls"], None),
    #cant test inductiva storage rm
    (["inductiva", "storage", "size"], None),
    (["inductiva", "task-runner"], None),
    #cant test inductiva task-runner launch
    #cant test inductiva task-runner remove
    # TODO: try to test inductiva task-runner on the github machines
    (["inductiva", "tasks"], None),
    (["inductiva", "tasks", "list"], None),
    #cant test inductiva tasks download
    #cant test inductiva tasks info
    #cant test inductiva tasks kill
    (["inductiva", "user"], None),
    (["inductiva", "user", "info"], None),
]


@pytest.mark.parametrize("command, user_input", CLI_COMMANDS)
def test_cli(command, user_input):
    process = subprocess.Popen(command,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               text=True)

    # Send user input if needed

    stdout, stderr = process.communicate(input=user_input, timeout=120)
    print("STDOUT: ", stdout)
    print("STDERR: ", stderr)
    print("#####################")

    assert process.returncode == 0
