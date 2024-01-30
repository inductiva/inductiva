import inductiva
import pytest
from unittest.mock import Mock
from inductiva import constants
from inductiva.client import exceptions
from inductiva.client.model.task_status_code import TaskStatusCode

from inductiva.client.paths.tasks_task_id_input.post import ApiResponseFor200


def test_task_kill__string_timeout__typeError_exception():
    """ 
    Check if task.kill returns TypeError when calling
    kill with timeout as a string.
    """
    task = inductiva.tasks.Task("123")

    with pytest.raises(TypeError) as exception:
        task.kill(wait_timeout="Olá mundo.")
    assert str(exception.value) == "Wait timeout must be a number."


def test_task_kill__negative_timeout__ValueError_exception():
    """ 
    Check if task.kill returns ValueError when calling
    kill with timeout as a negative number.
    """
    task = inductiva.tasks.Task("123")

    with pytest.raises(ValueError) as exception:
        task.kill(wait_timeout=-1)
    assert str(exception.value) == "Wait timeout must be a positive number."


def test_task_kill__zero_timeout__ValueError_exception():
    """
    Check if task.kill returns ValueError when calling
    kill with timeout as zero.
    """
    task = inductiva.tasks.Task("123")
    with pytest.raises(ValueError) as exception:
        task.kill(wait_timeout=0)
    assert str(exception.value) == "Wait timeout must be a positive number."


def test_task_kill__none_timeout__none():
    """
    Check if task.kill returns none when wait
    timeout is none.
    """
    task = inductiva.tasks.Task("123")
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
    task._send_kill_request = Mock(return_value=ApiResponseFor200)
    task._check_if_pending_kill = Mock(return_value=(pending_kill_success,
                                                     pending_kill_status))
    success = task.kill(wait_timeout=1)

    assert success == expected_output


@pytest.mark.parametrize("get_status_response", [
    (TaskStatusCode.PENDINGKILL),
    (TaskStatusCode.KILLED)(TaskStatusCode.ZOMBIE),
    (TaskStatusCode.PENDINGINPUT),
    (TaskStatusCode.FAILED),
])
def test__send_kill_request__positive_max_api_requests__none(
        get_status_response):
    """
    Check the _send_kill_request return for multiple
    get_status returns.
    """
    task = inductiva.tasks.Task("123")

    task._api.kill_task = Mock(return_value=ApiResponseFor200)
    task.get_status = Mock(return_value=get_status_response)

    kill_request_return = task._send_kill_request(
        constants.TASK_KILL_MAX_API_REQUESTS)

    assert kill_request_return is None


def test__send_kill_request__positive_max_api_requests__api_exception__runtimeError(
):
    """
    Check the _send_kill_request return when the api call
    gives an exception.
    """
    task = inductiva.tasks.Task("123")

    task._api.kill_task = Mock(
        side_effect=exceptions.ApiException(400, "Bad Request"))
    task.get_status = Mock(return_value=TaskStatusCode.PENDINGKILL)

    with pytest.raises(RuntimeError) as exception:
        _ = task._send_kill_request(constants.TASK_KILL_MAX_API_REQUESTS)
    assert str(
        exception.value
    ) == "Something went wrong while sending the kill command. Please try again later."


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
    task = inductiva.tasks.Task("123")
    task.get_status = Mock(return_value=status_code)
    success, status = task._check_if_pending_kill(2)

    assert success == expected_success and status == status_code
