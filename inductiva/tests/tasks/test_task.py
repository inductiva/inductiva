"""Test file for Tasks class."""
import inductiva
import pytest
from unittest.mock import Mock
from inductiva import constants
from inductiva.client import exceptions
from inductiva.client.model.task_status_code import TaskStatusCode

from inductiva.client.paths.tasks_task_id_input.put import ApiResponseFor200

import inductiva.client.paths.tasks_task_id_output_list.get as output_list_get
from inductiva.tasks.task import TaskInfo


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
    task = inductiva.tasks.Task("123")

    # pylint: disable=W0212
    task._api.kill_task = Mock(return_value=ApiResponseFor200)
    task.get_status = Mock(return_value=get_status_response)

    # pylint: disable=W0212
    kill_request_return = task._send_kill_request(
        constants.TASK_KILL_MAX_API_REQUESTS)

    assert kill_request_return is None


def test__send_kill_request__api_exception__runtimeerror():
    """
    Check the _send_kill_request return when the api call
    gives an exception.
    """
    task = inductiva.tasks.Task("123")

    # pylint: disable=W0212
    task._api.kill_task = Mock(
        side_effect=exceptions.ApiException(400, "Bad Request"))
    task.get_status = Mock(return_value=TaskStatusCode.PENDINGKILL)

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
    task = inductiva.tasks.Task("123")
    task.get_status = Mock(return_value=status_code)
    # pylint: disable=W0212
    success, status = task._check_if_pending_kill(2)

    assert success == expected_success and status == status_code


def test__get_output_info():
    """
    Check if the output info is correctly returned.
    """
    task = inductiva.tasks.Task("123")
    # pylint: disable=W0212
    mock_resp_body = output_list_get.SchemaFor200ResponseBodyApplicationJson(
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

    task._api.get_outputs_list = Mock(
        return_value=output_list_get.ApiResponseFor200(
            response=Mock(),
            body=mock_resp_body,
        ))

    output_info = task.get_output_info()

    assert output_info.n_files == 2
    assert output_info.total_size_bytes == 320
    assert output_info.total_compressed_size_bytes == 150


task_info_dic = {
    "task_id": "gfuhdirkeqypp3b5cb6ped1re",
    "status": "success",
    "method_name": "amrWind.amrWind.run_simulation",
    "storage_path": None,
    "container_image": "docker://inductiva/kutu:amr-wind_v1.4.0",
    "project": "pbarbosaaa788d1b",
    "create_time": "2024-07-12T09:57:58.111714+00:00",
    "input_submit_time": "2024-07-12T09:57:58.664646+00:00",
    "start_time": "2024-07-12T09:57:58.666103+00:00",
    "computation_start_time": "2024-07-12T09:57:58.773196+00:00",
    "computation_end_time": "2024-07-12T09:58:03.895437+00:00",
    "end_time": "2024-07-12T09:58:05.868836+00:00",
    "cost": 0.00016,
    "storage_size": 12763590,
    "metrics": {
        "total_seconds": 7.757,
        "container_image_download_seconds": None,
        "queue_time_seconds": 0.001,
        "computation_seconds": 5.122,
        "input_upload_seconds": 0.55,
        "input_download_seconds": 0.081,
        "input_decompression_seconds": 0.003,
        "output_compression_seconds": 1.63,
        "output_upload_seconds": 0.309,
        "input_zipped_size_bytes": 1437,
        "input_size_bytes": 11294,
        "output_total_files": 92,
        "output_size_bytes": 52241441,
        "output_zipped_size_bytes": 12762153
    },
    "executer": {
        "uuid": "d53a581f-f796-4e8d-8e93-104fc747352a",
        "cpu_count_logical": 4,
        "cpu_count_physical": 2,
        "memory": 16785870848,
        "n_mpi_hosts": 0,
        "vm_type": "e2-standard-4",
        "vm_name": "default-machine-group-dm1l",
        "host_type": "GCP",
        "error_detail": None
    }
}


def test_taskinfo_innit():

    task_info = TaskInfo(**task_info_dic)

    assert task_info.task_id == "gfuhdirkeqypp3b5cb6ped1re"


def test_taskinfo_format_time_metric():

    task_info = TaskInfo(**task_info_dic)
    # pylint: disable=W0212
    result_smaller_60 = task_info._format_time_metric(
        "total_seconds", task_info.time_metrics.total_seconds.value)

    result_bigger_60 = task_info._format_time_metric("total_seconds", 61.5)

    task_info.is_running = True

    result_still_running = task_info._format_time_metric(
        "computation_seconds", 61.5)

    task_info.is_terminal = True
    result_container = task_info._format_time_metric(
        "container_image_download_seconds", None)
    result_comp_seconds = task_info._format_time_metric("computation_seconds",
                                                        None)

    task_info.is_terminal = False
    result_container_na = task_info._format_time_metric(
        "container_image_download_seconds", None)

    result_compression = task_info._format_time_metric(
        "output_compression_seconds", None)

    result_final_return = task_info._format_time_metric("hello_world", None)

    assert (result_smaller_60 == "7.76 s" and \
            result_bigger_60 == "0:01:02" and \
            "still running" in result_still_running and \
            "used cached" in result_container and \
            "N/A" in result_comp_seconds and \
            "N/A" in result_container_na and \
            "until task ends" in result_compression and \
            result_final_return is None
            )
