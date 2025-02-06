"""Test for the CLI commands."""
import subprocess
import pytest
import os


#read INDUCTIVA_API_KEY from os env vars
API_KEY = os.environ["INDUCTIVA_API_KEY"]

# ([command], user_input, can_timeout)
CLI_COMMANDS = [
    (["inductiva"], None, False),
    (["inductiva", "auth"], None, False),
    (["inductiva", "auth", "login"], API_KEY, False),
    (["inductiva", "auth", "logout"], None, False),

    #cant test inductiva logs
    #(["inductiva", "logs"], None, True),
    (["inductiva", "projects"], None, False),
    (["inductiva", "projects", "ls"], None, False),
    (["inductiva", "quotas"], None, False),
    (["inductiva", "quotas", "list"], None, False),
    (["inductiva", "resources"], None, False),
    ([
        "inductiva",
        "resources",
        "start",
        "c2-standard-4",
    ], None, False),
    (["inductiva", "resources", "cost", "c2-standard-4"], None, False),
    #cant test inductiva resources info
    (["inductiva", "resources", "available"], None, False),
    (["inductiva", "resources", "ls"], None, False),
    (["inductiva", "resources", "terminate", "--all"], "yes", False),
    (["inductiva", "simulators"], None, False),
    (["inductiva", "simulators", "list"], None, False),
    (["inductiva", "storage"], None, False),
    (["inductiva", "storage", "ls"], None, False),
    #cant test inductiva storage rm
    (["inductiva", "storage", "size"], None, False),
    (["inductiva", "task-runner"], None, False),
    #cant test inductiva task-runner launch
    #cant test inductiva task-runner remove
    # TODO: try to test inductiva task-runner on the github machines
    (["inductiva", "tasks"], None, False),
    (["inductiva", "tasks", "list"], None, False),
    #cant test inductiva tasks download
    #cant test inductiva tasks info
    #cant test inductiva tasks kill
    (["inductiva", "user"], None, False),
    (["inductiva", "user", "info"], None, False),
]


@pytest.mark.parametrize("command, user_input, can_timeout", CLI_COMMANDS)
def test_cli(command, user_input, can_timeout):
    process = subprocess.Popen(command,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               text=True)

    # Send user input if needed
    try:
        stdout, stderr = process.communicate(input=user_input, timeout=120)
        print("STDOUT: ", stdout)
        print("STDERR: ", stderr)
        print("#####################")
    except subprocess.TimeoutExpired:
        if can_timeout:
            process.returncode = 0
        print("TimeoutExpired")
        print("#####################")

    assert process.returncode == 0
