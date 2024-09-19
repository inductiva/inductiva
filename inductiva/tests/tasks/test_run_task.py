"""Test file for the run_simulation file."""
import uuid
import os
import json
import tempfile

import pytest
from unittest import mock

import inductiva
from inductiva import tasks

# pylint: disable=R1732
TEMP_DIR = tempfile.TemporaryDirectory()
TASK_METADATA_FILENAME = "task_metadata.json"


# pylint: disable=unused-argument
def _api_invoker(task_id, *args, **kwargs):
    """Dummy api_invoker argument for testing. Returns the api_method_name
    """
    return task_id


def _id_in_metadata_file(task_id):
    """Checks if the ID exists in the file."""
    if not os.path.exists(TASK_METADATA_FILENAME):
        return False
    with open(TASK_METADATA_FILENAME, "r", encoding="utf-8") as f:
        json_list = [json.loads(line) for line in f]
    return task_id in map(lambda x: x["task_id"], json_list)


@pytest.mark.parametrize("task_id,disable_logging", [(".id_1", True),
                                                     (".id_2", False),
                                                     (".id_3", True),
                                                     (".id_4", False)])
def test_run_simulation_logging(task_id, disable_logging):
    """Tests if the id of the task was added to the file."""
    inductiva.set_output_dir(None)
    os.chdir(TEMP_DIR.name)
    dummy_dir = os.path.join(os.getcwd(), f"dummy_dir_{task_id}")
    if not os.path.exists(dummy_dir):
        os.mkdir(dummy_dir)
    os.environ["DISABLE_TASK_METADATA_LOGGING"] = str(disable_logging)
    inductiva.set_api_key("DUMMY")

    mock_mg = mock.Mock()
    mock_mg.id = str(uuid.uuid4())

    tasks.run_simulation(api_method_name=task_id,
                         input_dir=dummy_dir,
                         computational_resources=mock_mg,
                         api_invoker=_api_invoker)
    assert _id_in_metadata_file(task_id) == (not disable_logging)
