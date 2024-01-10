"""Test file for the run_simulation file."""
import os
import json

import pathlib

from inductiva import tasks

LOGGED_TASK_ID = "1234"
NOT_LOGGED_TASK_ID = "4321"


def _id_in_metadata_file(task_id):
    """Checks if the ID exists in the file."""
    with open(tasks.TASK_METADATA_FILENAME, "r", encoding="utf-8") as f:
        json_list = [json.loads(line) for line in f]
    return task_id in map(lambda x: x["task_id"], json_list)


def test_run_simulation_logging(mocker):
    """Tests if the id of the task was added to the file."""

    # Old logs may cause problems so we remove the file.
    if os.path.exists(tasks.TASK_METADATA_FILENAME):
        os.remove(tasks.TASK_METADATA_FILENAME)

    invoke_api_mock = mocker.patch("inductiva.api.methods.invoke_async_api")
    mocker.patch("inductiva.tasks.Task")

    input_dir = pathlib.Path().cwd()

    _run_simulation_and_assert("some_api_method",
                               input_dir,
                               invoke_api_mock,
                               LOGGED_TASK_ID,
                               disable_logging=False,
                               id_exists_in_file=True)
    _run_simulation_and_assert("some_api_method",
                               input_dir,
                               invoke_api_mock,
                               NOT_LOGGED_TASK_ID,
                               disable_logging=True,
                               id_exists_in_file=False)


def _run_simulation_and_assert(api_method_name, input_dir, invoke_api_mock,
                               task_id, disable_logging, id_exists_in_file):
    os.environ["DISABLE_TASK_METADATA_LOGGING"] = str(disable_logging)
    invoke_api_mock.return_value = task_id
    tasks.run_simulation(api_method_name=api_method_name, input_dir=input_dir)
    assert _id_in_metadata_file(task_id) == id_exists_in_file
