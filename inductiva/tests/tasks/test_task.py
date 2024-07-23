"""Test file for Tasks class."""
import inductiva
import pytest
from unittest.mock import Mock, patch
from inductiva import constants
from inductiva.client import exceptions
import inductiva.client
from inductiva.client.model.task_status_code import TaskStatusCode

import inductiva.client.paths
from inductiva.client.paths.tasks_task_id_input.put import ApiResponseFor200

import inductiva.client.paths.tasks_task_id_kill
import inductiva.client.paths.tasks_task_id_kill.post
import inductiva.client.paths.tasks_task_id_output_list.get as output_list_get


def test_task_kill__string_timeout__typeerror_exception():
    """
    Check if task.kill returns TypeError when calling
    kill with timeout as a string.
    """
    task = inductiva.tasks.Task("123")

    with pytest.raises(TypeError) as exception:
        task.kill(wait_timeout="Olá mundo.")
    assert "number" in str(exception.value)


def test_task_kill__negative_timeout__valueerror_exception():
    """
    Check if task.kill returns ValueError when calling
    kill with timeout as a negative number.
    """
    task = inductiva.tasks.Task("123")

    with pytest.raises(ValueError) as exception:
        task.kill(wait_timeout=-1)
    assert "positive number" in str(exception.value)


def test_task_kill__zero_timeout__valueerror_exception():
    """
    Check if task.kill returns ValueError when calling
    kill with timeout as zero.
    """
    task = inductiva.tasks.Task("123")
    with pytest.raises(ValueError) as exception:
        task.kill(wait_timeout=0)
    assert "positive number" in str(exception.value)


def test_task_kill__none_timeout__none():
    """
    Check if task.kill returns none when wait
    timeout is none.
    """
    task = inductiva.tasks.Task("123")
    # pylint: disable=W0212
    task._send_kill_request = Mock(return_value=None)
    kill_return = task.kill()

    assert kill_return is None


@pytest.mark.parametrize(
    "pending_kill_success, pending_kill_status, expected_output", [
        (True, TaskStatusCode.KILLED, True),
        (False, TaskStatusCode.PENDINGKILL, False),
        (True, TaskStatusCode.ZOMBIE, False),
        (True, TaskStatusCode.FAILED, False),
        (True, TaskStatusCode.PENDINGINPUT, False),
    ])
def test_task_kill__positive_timeout__success(pending_kill_success,
                                              pending_kill_status,
                                              expected_output):
    """
    Check task.kill return based on the output of
    _check_if_pending_kill.
    """
    task = inductiva.tasks.Task("123")
    # pylint: disable=W0212
    task._send_kill_request = Mock(return_value=ApiResponseFor200)
    # pylint: disable=W0212
    task._check_if_pending_kill = Mock(return_value=(pending_kill_success,
                                                     pending_kill_status))
    success = task.kill(wait_timeout=1)

    assert success == expected_output


@pytest.mark.parametrize("get_status_response", [
    (TaskStatusCode.PENDINGKILL),
    (TaskStatusCode.KILLED),
    (TaskStatusCode.ZOMBIE),
    (TaskStatusCode.PENDINGINPUT),
    (TaskStatusCode.FAILED),
])
def test__send_kill_request__positive_max_api_requests__none(
        get_status_response):
    """
    Check the _send_kill_request return for multiple
    get_status returns.
    """

    # Some Mock classes
    class TaskApiMock:
        """Some Mock Class"""

        def __init__(self, *_args, **_kwargs):
            pass

        def kill_task(self, *_args, **_kwargs):
            return

    with patch("inductiva.tasks.task.tasks_api") as mock_api, patch(
            "inductiva.tasks.Task.get_status",
            return_value=get_status_response), patch(
                "inductiva.tasks.Task.is_terminal", return_value=False):
        mock_api.TasksApi = TaskApiMock

        task = inductiva.tasks.Task("123")
        # pylint: disable=W0212
        kill_request_return = task._send_kill_request(
            constants.TASK_KILL_MAX_API_REQUESTS)

        assert kill_request_return is None


def test__send_kill_request__api_exception__runtimeerror():
    """
    Check the _send_kill_request return when the api call
    gives an exception.
    """

    # Some Mock classes
    class TaskApiMock:
        """Some Mock Class"""

        def __init__(self, *_args, **_kwargs):
            pass

        def kill_task(self, *_args, **_kwargs):
            raise exceptions.ApiException(400, "Bad Request")

    with patch("inductiva.tasks.task.tasks_api") as mock_api, patch(
            "inductiva.tasks.Task.get_status",
            return_value=TaskStatusCode.PENDINGKILL), patch(
                "inductiva.tasks.Task.is_terminal", return_value=False), patch(
                    "inductiva.constants.TASK_KILL_RETRY_SLEEP_SEC", 0.1):
        mock_api.TasksApi = TaskApiMock

        task = inductiva.tasks.Task("123")

        with pytest.raises(RuntimeError) as exception:
            # pylint: disable=W0212
            _ = task._send_kill_request(constants.TASK_KILL_MAX_API_REQUESTS)
        assert "kill command" in str(exception.value)


@pytest.mark.parametrize("status_code, expected_success", [
    (TaskStatusCode.PENDINGKILL, False),
    (TaskStatusCode.KILLED, True),
    (TaskStatusCode.ZOMBIE, True),
    (TaskStatusCode.FAILED, True),
])
def test__check_if_pending_kill__wait_timeout_positive__success_status(
        status_code, expected_success):
    """
    Check the _check_if_pending_kill returns for the different satus
    returned by the api.
    """
    with patch("inductiva.tasks.Task.get_status",
               return_value=status_code), patch(
                   "inductiva.constants.TASK_KILL_RETRY_SLEEP_SEC", 0.1):
        task = inductiva.tasks.Task("123")
        # pylint: disable=W0212

        success, status = task._check_if_pending_kill(0.1)

        assert success == expected_success and status == status_code


def test__get_output_info():
    """
    Check if the output info is correctly returned.
    """

    # Some Mock classes
    class TaskApiMock:
        """Some Mock Class"""

        def __init__(self, *_args, **_kwargs):
            pass

        def get_outputs_list(self, *_args, **_kwargs):
            body = output_list_get.SchemaFor200ResponseBodyApplicationJson(
                size=320,
                contents=[
                    {
                        "name": "file1.txt",
                        "size": 100,
                        "compressed_size": 50
                    },
                    {
                        "name": "file2.txt",
                        "size": 200,
                        "compressed_size": 100
                    },
                ],
            )
            return output_list_get.ApiResponseFor200(
                response=Mock(),
                body=body,
            )

    with patch("inductiva.tasks.task.tasks_api") as mock_api:
        mock_api.TasksApi = TaskApiMock

        task = inductiva.tasks.Task("123")
        output_info = task.get_output_info()

        assert output_info.n_files == 2
        assert output_info.total_size_bytes == 320
        assert output_info.total_compressed_size_bytes == 150
