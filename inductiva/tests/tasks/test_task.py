"""Test file for Tasks class."""
from unittest import mock

import inductiva
import pytest
from unittest.mock import Mock
from inductiva import constants
from inductiva.client import exceptions
import inductiva.client
from inductiva.client.model.task_status_code import TaskStatusCode

import inductiva.client.paths
import inductiva.client.paths.tasks_task_id
from inductiva.client.paths.tasks_task_id_input.put import ApiResponseFor200

import inductiva.client.paths.tasks_task_id_output_list.get as output_list_get
import inductiva.client.paths.tasks_task_id_position_in_queue
import inductiva.client.paths.tasks_task_id_position_in_queue.get
import inductiva.client.paths.tasks_task_id_status
import inductiva.client.paths.tasks_task_id_status.get
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


def test_task_kill__wrong_verbosity_level():
    """
    Check if task.kill throws an exception when verbosity level is not valid.
    """
    task = inductiva.tasks.Task("123")

    with pytest.raises(ValueError) as exception:
        task.kill(verbosity_level=3)
    assert "level not allowed" in str(exception.value)


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
        "error_detail": "test error detail"
    }
}

task_info_bytes = (
    b'{"task_id":"6gpbvrxr46dvm8p7i1hcjt419",'
    b'"status":"success","method_name":"amrWind.amrWind.run_simulation",'
    b'"storage_path":null,"container_image":'
    b'"docker://inductiva/kutu:amr-wind_v1.4.0","project":"user1",'
    b'"create_time":"2024-07-16T15:07:58.211579+00:00","input_submit_time"'
    b':"2024-07-16T15:07:58.817154+00:00","start_time":'
    b'"2024-07-16T15:07:58.819051+00:00","computation_start_time":'
    b'"2024-07-16T15:07:58.998724+00:00","computation_end_time":'
    b'"2024-07-16T15:08:05.166285+00:00","end_time":'
    b'"2024-07-16T15:08:06.932283+00:00","cost":0.0082,"storage_size":1243,'
    b'"metrics":{"total_seconds":8,"container_image_download_seconds":null,'
    b'"queue_time_seconds":0.002,"computation_seconds":6.168,'
    b'"input_upload_seconds":0.603,"input_download_seconds":0.084,'
    b'"input_decompression_seconds":0.001,"output_compression_seconds":1.4,'
    b'"output_upload_seconds":0.248,"input_zipped_size_bytes":1437,'
    b'"input_size_bytes":11294,"output_total_files":44,"output_size_bytes":'
    b'52241740,"output_zipped_size_bytes":12747906},"executer":{"uuid":'
    b'"4dbc0ffc-5d23-4365-88b9-7571be845dd3","cpu_count_logical":32,'
    b'"cpu_count_physical":16,"memory":135081934,"n_mpi_hosts":0,"vm_type":'
    b'"c2d-standard-32","vm_name":"api-tfmpeknbb41jc6k6n31haewff-7vz2",'
    b'"host_type":"GCP","error_detail":null}}')


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


def test_format_data_metric():
    task_info = TaskInfo(**task_info_dic)
    # pylint: disable=W0212
    result_bytes = task_info._format_data_metric("input_size_bytes", 123)
    result_not_bytes = task_info._format_data_metric("hello world", 123)

    assert result_bytes == "123 B" and result_not_bytes == 123


def test_repr_():
    task_info = TaskInfo(**task_info_dic)
    # pylint: disable=c2801
    result = task_info.__repr__()
    assert ("Container image download:  N/A (used cached image)" in result and
            "52.24 MB" in result and "test error detail" in result)


def test_str_():
    task_info = TaskInfo(**task_info_dic)
    # pylint: disable=c2801
    result = task_info.__str__()
    print(result)
    assert ("Container image download:  N/A (used cached image)" in result and
            "52.24 MB" in result and "test error detail" in result)


@pytest.mark.parametrize("get_status_response,expected_result", [
    (TaskStatusCode.PENDINGINPUT, False),
    (TaskStatusCode.PENDINGKILL, False),
    (TaskStatusCode.STARTED, True),
    (TaskStatusCode.KILLED, False),
    (TaskStatusCode.ZOMBIE, False),
    (TaskStatusCode.FAILED, False),
])
def test_is_running(get_status_response, expected_result):
    task = inductiva.tasks.Task("123")
    with mock.patch("inductiva.tasks.Task.get_status",
                    return_value=get_status_response):
        result = task.is_running()
        assert result == expected_result


@pytest.mark.parametrize("get_status_response,expected_result", [
    (TaskStatusCode.PENDINGINPUT, False),
    (TaskStatusCode.PENDINGKILL, False),
    (TaskStatusCode.STARTED, False),
    (TaskStatusCode.KILLED, True),
    (TaskStatusCode.ZOMBIE, True),
    (TaskStatusCode.FAILED, True),
])
def test_is_failed(get_status_response, expected_result):
    task = inductiva.tasks.Task("123")
    with mock.patch("inductiva.tasks.Task.get_status",
                    return_value=get_status_response):
        result = task.is_failed()
        assert result == expected_result


def test_from_api_info():
    task = inductiva.tasks.Task.from_api_info(task_info_dic)
    # pylint: disable=W0212
    assert (task.id == task_info_dic["task_id"] and
            task._status == task_info_dic["status"] and
            task.info.start_time == task_info_dic["start_time"])


@pytest.mark.parametrize("status, status_code, tasks_ahead", [
    ("pending-input", TaskStatusCode.PENDINGINPUT, 1),
    ("pending-kill", TaskStatusCode.PENDINGKILL, 2),
    ("started", TaskStatusCode.STARTED, 3),
    ("killed", TaskStatusCode.KILLED, 4),
    ("zombie", TaskStatusCode.ZOMBIE, 5),
    ("failed", TaskStatusCode.FAILED, 0),
])
def test_get_status(status, status_code, tasks_ahead):
    get_task_status_return = mock.MagicMock()
    get_task_status_return.body = {
        "id": "6gpbvrxr46dvm8p7i1hcjt419",
        "status": status,
        "position_in_queue": {
            "tasks_ahead": tasks_ahead
        }
    }
    with mock.patch(
            "inductiva.client.paths.tasks_task_id_status.get."
            "GetTaskStatus.get_task_status",
            return_value=get_task_status_return):
        task = inductiva.tasks.Task("123")
        result = task.get_status()
        # pylint: disable=W0212
        assert (task._tasks_ahead == tasks_ahead and result == status_code)


def test_get_status_terminal():
    task = inductiva.tasks.Task("123")
    # pylint: disable=W0212
    task._status = TaskStatusCode.KILLED
    result = task.get_status()
    assert result == TaskStatusCode.KILLED


@pytest.mark.parametrize("tasks_ahead", [
    (1),
    (2),
    (3),
    (4),
    (0),
])
def test_get_position_in_queue(tasks_ahead):
    get_task_position_return = mock.MagicMock()
    get_task_position_return.body = {"tasks_ahead": tasks_ahead}
    with mock.patch(
            "inductiva.client.paths.tasks_task_id_position_in_queue.get"
            ".GetTaskPositionInQueue.get_task_position_in_queue",
            return_value=get_task_position_return):
        task = inductiva.tasks.Task("123")
        result = task.get_position_in_queue()
        assert result == tasks_ahead


def test_get_position_in_queue_exception():
    with mock.patch(
            "inductiva.client.paths.tasks_task_id_position_in_queue.get"
            ".GetTaskPositionInQueue.get_task_position_in_queue",
            side_effect=exceptions.ApiException(404, "Not Found")):
        task = inductiva.tasks.Task("123")
        result = task.get_position_in_queue()
        assert result is None


def test_get_info():
    get_task_return = mock.MagicMock()
    get_task_return.response.data = task_info_bytes
    with mock.patch("inductiva.client.paths.tasks_task_id.get.GetTask.get_task",
                    return_value=get_task_return):
        task = inductiva.tasks.Task("123")
        info = task.get_info()

        assert info.task_id == "6gpbvrxr46dvm8p7i1hcjt419"


def test_info_property_none():
    get_task_return = mock.MagicMock()
    get_task_return.response.data = task_info_bytes
    with mock.patch("inductiva.client.paths.tasks_task_id.get.GetTask.get_task",
                    return_value=get_task_return):
        task = inductiva.tasks.Task("123")
        # pylint: disable=W0212
        task._info = None
        info = task.info
        assert info.task_id == "6gpbvrxr46dvm8p7i1hcjt419"


@pytest.mark.parametrize("tasks_ahead, is_tty, expected_output", [
    (0, True, "The task 123 is about to start."),
    (0, False,
     "The task 123 is about to start.                                          "
    ), (1, True, "Number of tasks ahead of task 123 in queue: 1"),
    (1, False,
     "Number of tasks ahead of task 123 in queue: 1                            "
    )
])
def test_setup_queue_message(tasks_ahead, is_tty, expected_output):
    task = inductiva.tasks.Task("123")
    # pylint: disable=W0212
    task._tasks_ahead = tasks_ahead
    result = task._setup_queue_message(is_tty=is_tty)
    print(result)
    print(expected_output)
    assert result == expected_output


#Wait is realy hard to test since it depends on the status changing to avoid
#an infinite loop. We will test only the cases where the status is terminal
@pytest.mark.parametrize(
    "get_status_response",
    [
        (TaskStatusCode.KILLED),
        (TaskStatusCode.ZOMBIE),
        (TaskStatusCode.FAILED),
        #This status will make an api call and we already have a test for that
        #(TaskStatusCode.EXECUTERFAILED),
        (TaskStatusCode.EXECUTERTERMINATED),
        (TaskStatusCode.EXECUTERTERMINATEDBYUSER),
        (TaskStatusCode.SPOTINSTANCEPREEMPTED),
        (TaskStatusCode.EXECUTERTERMINATEDTTLEXCEEDED),
        (TaskStatusCode.TTLEXCEEDED),
        (TaskStatusCode.SUCCESS),
        #All non terminal statuses will result in a infinite loop inside the wait
    ])
def test_wait(get_status_response):
    task = inductiva.tasks.Task("123")
    with mock.patch("inductiva.tasks.Task.get_status",
                    return_value=get_status_response):
        result = task.wait()
        assert result == get_status_response
