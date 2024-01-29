"""Tests on a dummy simulator for testing purposes."""
import tempfile
import json
import os

from inductiva.simulators.dummy_simulator import DummySimulator

INPUT_DIR = os.path.join(os.path.dirname(__file__), "test_input_dir")
BACKEND_ARGS_FILE = "test_arguments.json"
TEMP_DIR = tempfile.TemporaryDirectory()


def validate_args_from_backend(test_args, output_dir):
    """Validate that the arguments are sent and validated correctly."""
    backend_args_file = os.path.join(output_dir, BACKEND_ARGS_FILE)
    with open(backend_args_file, "r", encoding="utf-8") as f:
        backend_args = json.load(f)

    assert test_args["input_filename"] == backend_args["input_filename"]
    assert set(test_args["input_dir_list"]).issubset(
        set(backend_args["input_dir_list"]))
    assert test_args["sleep_time"] == backend_args["sleep_time"]
    assert test_args["commands"] == backend_args["commands"]


def test_dummy_simulator__send_args__validate_args():
    """Test that the arguments are sent and validated correctly."""

    test_args = {
        "input_dir_list": os.listdir(INPUT_DIR),
        "input_filename": "who_am_i.txt",
        "commands": None,
        "sleep_time": 1
    }

    dummy = DummySimulator()

    task = dummy.run(input_dir=INPUT_DIR,
                     input_filename=test_args["input_filename"],
                     sleep_time=test_args["sleep_time"])

    task.wait()
    output_dir = task.download_outputs(output_dir=str(TEMP_DIR))

    assert os.path.exists(os.path.join(output_dir, BACKEND_ARGS_FILE))
    validate_args_from_backend(test_args, output_dir)


def test_dummy_simulator__with_commands():

    test_args = {
        "input_dir_list": os.listdir(INPUT_DIR),
        "input_filename": "who_am_i.txt",
        "commands": [{
            "cmd": "chmod +x logs.sh",
            "prompts": []
        }, {
            "cmd": "./logs.sh",
            "prompts": []
        }],
        "sleep_time": 1
    }

    # Initialize the Simulator
    dummy = DummySimulator()

    task = dummy.run(input_dir=INPUT_DIR,
                     input_filename=test_args["input_filename"],
                     commands=test_args["commands"],
                     sleep_time=test_args["sleep_time"])

    task.wait()
    output_dir = task.download_outputs(output_dir=str(TEMP_DIR))

    validate_args_from_backend(test_args, output_dir)

    assert os.path.exists(os.path.join(output_dir, "logs.txt"))
    assert os.path.exists(os.path.join(output_dir, "stdout.txt"))
    assert os.path.exists(os.path.join(output_dir, "stderr.txt"))
