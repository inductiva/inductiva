"""Tests for the benchmarks module."""
from unittest import mock

from inductiva.benchmarks.methods import _render_dict
from inductiva.benchmarks.methods import _compute_tasks_to_run
from inductiva.benchmarks.methods import _can_start_resource
from inductiva.benchmarks.methods import _tasks_by_vm_type
from inductiva.tasks.task import TaskInfo

executer_4 = mock.MagicMock()
executer_4.vm_type = "c2-standard-4"
executer_8 = mock.MagicMock()
executer_8.vm_type = "c2-standard-8"
executer_4d = mock.MagicMock()
executer_4d.vm_type = "c2d-standard-4"
executer_8d = mock.MagicMock()
executer_8d.vm_type = "c2d-standard-8"


def test___tasks_by_vm_type():
    """Test the tasks_by_vm_type function
    when a project has multiple vm_types.
    """
    task_info_4 = TaskInfo()
    task_info_4.executer = executer_4d
    task_info_8 = TaskInfo()
    task_info_8.executer = executer_8d

    project = mock.MagicMock()
    task1 = mock.MagicMock()
    task1.get_info = mock.MagicMock(return_value=task_info_4)
    task2 = mock.MagicMock()
    task2.get_info = mock.MagicMock(return_value=task_info_8)
    task3 = mock.MagicMock()
    task3.get_info = mock.MagicMock(return_value=task_info_4)

    project.get_tasks = mock.MagicMock(return_value=[task1, task2, task3])
    # pylint: disable=W0212
    current_project_tasks = _tasks_by_vm_type(project)

    assert current_project_tasks == {
        "c2d-standard-4": [task1, task3],
        "c2d-standard-8": [task2]
    }


def test___tasks_by_vm_type__no_tasks():
    """Test the tasks_by_vm_type function when a project has no tasks."""
    project = mock.MagicMock()

    project.get_tasks = mock.MagicMock(return_value=[])
    # pylint: disable=W0212
    current_project_tasks = _tasks_by_vm_type(project)

    assert not current_project_tasks


def test___tasks_by_vm_type__all_same_vm_type():
    """Test the tasks_by_vm_type function
    when all tasks have the same vm_type.
    """
    task_info_4 = TaskInfo()
    task_info_4.executer = executer_4d

    project = mock.MagicMock()
    task1 = mock.MagicMock()
    task1.get_info = mock.MagicMock(return_value=task_info_4)
    task2 = mock.MagicMock()
    task2.get_info = mock.MagicMock(return_value=task_info_4)
    task3 = mock.MagicMock()
    task3.get_info = mock.MagicMock(return_value=task_info_4)

    project.get_tasks = mock.MagicMock(return_value=[task1, task2, task3])
    # pylint: disable=W0212
    current_project_tasks = _tasks_by_vm_type(project)

    assert current_project_tasks == {"c2d-standard-4": [task1, task2, task3]}


def test__compute_tasks_to_run_some_task_ran():
    """Test the _compute_tasks_to_run function when some task already ran."""
    machines = ["c2-standard-4", "c2-standard-8"]

    task_info_4 = TaskInfo()
    task_info_4.executer = executer_4
    task_info_8 = TaskInfo()
    task_info_8.executer = executer_8

    task1 = mock.MagicMock()
    task1.get_info = mock.MagicMock(return_value=task_info_4)
    task2 = mock.MagicMock()
    task2.get_info = mock.MagicMock(return_value=task_info_8)
    task3 = mock.MagicMock()
    task3.get_info = mock.MagicMock(return_value=task_info_4)
    current_project_tasks = {
        "c2-standard-4": [task1, task3],
        "c2-standard-8": [task2]
    }
    replicas = 6
    # pylint: disable=W0212
    to_run = _compute_tasks_to_run(machines, current_project_tasks, replicas)
    assert to_run == {"c2-standard-4": 4, "c2-standard-8": 5}


def test__compute_tasks_to_run_no_task_ran():
    """Test the _compute_tasks_to_run function when no task have ran."""
    machines = ["c2-standard-4", "c2-standard-8"]
    current_project_tasks = {}
    replicas = 6
    # pylint: disable=W0212
    to_run = _compute_tasks_to_run(machines, current_project_tasks, replicas)
    assert to_run == {"c2-standard-4": 6, "c2-standard-8": 6}


def test__compute_tasks_to_run_all_task_ran():
    """Test the _compute_tasks_to_run function when all tasks have ran."""
    machines = ["c2-standard-4", "c2-standard-8"]

    task_info_4 = TaskInfo()
    task_info_4.executer = executer_4

    task_info_8 = TaskInfo()
    task_info_8.executer = executer_8

    task1 = mock.MagicMock()
    task1.get_info = mock.MagicMock(return_value=task_info_4)
    task2 = mock.MagicMock()
    task2.get_info = mock.MagicMock(return_value=task_info_8)
    task3 = mock.MagicMock()
    task3.get_info = mock.MagicMock(return_value=task_info_4)
    task4 = mock.MagicMock()
    task4.get_info = mock.MagicMock(return_value=task_info_8)
    current_project_tasks = {
        "c2-standard-4": [task1, task3],
        "c2-standard-8": [task2, task4]
    }
    replicas = 2
    # pylint: disable=W0212
    to_run = _compute_tasks_to_run(machines, current_project_tasks, replicas)

    assert {'c2-standard-4': 0, 'c2-standard-8': 0} == to_run


def test__render_dict_all_callable():
    """Test the _render_dict function
    when all values are callable.
    """
    dictionary = {
        "a": lambda x: x + 1,
        "b": lambda x: x + 2,
        "c": lambda x: x + 3
    }
    # pylint: disable=W0212
    result = _render_dict(dictionary, 1)
    assert result == {"a": 2, "b": 3, "c": 4}


def test__render_dict_no_callable():
    """Test the _render_dict function
    when no values are callable.
    """
    dictionary = {"a": 1, "b": 1, "c": 3}
    # pylint: disable=W0212
    result = _render_dict(dictionary, 1)
    assert result == {"a": 1, "b": 1, "c": 3}


def test__render_dict_some_callable():
    """Test the _render_dict function
    when some values are callable.
    """
    dictionary = {"a": 1, "b": lambda x: x + 20, "c": 3}
    # pylint: disable=W0212
    result = _render_dict(dictionary, 1)
    assert result == {"a": 1, "b": 21, "c": 3}


def test__render_dict_argument_is_object():
    """Test the _render_dict function
    when the argument passed is an object.
    """
    dictionary = {"a": 1, "b": lambda x: x.hello(), "c": 3}

    argument = mock.MagicMock()
    argument.hello = mock.MagicMock(return_value="hello world")

    # pylint: disable=W0212
    result = _render_dict(dictionary, argument)
    assert result == {"a": 1, "b": "hello world", "c": 3}


def fake_get_quotas_small():
    """Fake get_quotas function."""
    return {
        "max_vcpus": {
            "max_allowed": 1.0,
            "in_use": 0.0
        },
        "max_instances": {
            "max_allowed": 1.0,
            "in_use": 0.0
        },
        "max_price_hour": {
            "max_allowed": 1.0,
            "in_use": 0.0
        },
        "max_disk_size": {
            "max_allowed": 1.0,
            "in_use": None
        },
        "mg_max_idle": {
            "max_allowed": 1.0,
            "in_use": None
        },
        "mg_max_ttl": {
            "max_allowed": 1.0,
            "in_use": None
        }
    }


def fake_get_quotas_normal():
    """Fake get_quotas function."""
    return {
        "max_vcpus": {
            "max_allowed": 100.0,
            "in_use": 0.0
        },
        "max_instances": {
            "max_allowed": 10.0,
            "in_use": 0.0
        },
        "max_price_hour": {
            "max_allowed": 20.0,
            "in_use": 0.0
        },
        "max_disk_size": {
            "max_allowed": 50.0,
            "in_use": None
        },
        "mg_max_idle": {
            "max_allowed": 150.0,
            "in_use": None
        },
        "mg_max_ttl": {
            "max_allowed": 48.0,
            "in_use": None
        }
    }


def fake_get_quotas_no_total_num_machines():
    """Fake get_quotas function."""
    return {
        "max_vcpus": {
            "max_allowed": 100.0,
            "in_use": 0.0
        },
        "max_instances": {
            "max_allowed": 0.0,
            "in_use": 0.0
        },
        "max_price_hour": {
            "max_allowed": 20.0,
            "in_use": 0.0
        },
        "max_disk_size": {
            "max_allowed": 50.0,
            "in_use": None
        },
        "mg_max_idle": {
            "max_allowed": 150.0,
            "in_use": None
        },
        "mg_max_ttl": {
            "max_allowed": 48.0,
            "in_use": None
        }
    }


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_small)
def test__can_start_resource_cant_start_low_quotas():
    """Test the _can_start_resource function when the quotas are small."""
    resource = mock.MagicMock()
    resource.machine_type = "c2-standard-4"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.04724)
    resource.consumed_vcpus = 4

    res = _can_start_resource(resource)

    assert not res


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_normal)
def test__can_start_resource_cant_start_too_many_vcpu():
    """Test the _can_start_resource function when the resource asks for
    too many vcpus."""

    resource = mock.MagicMock()
    resource.machine_type = "c2-standard-400"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.04724)
    resource.consumed_vcpus = 400

    res = _can_start_resource(resource)

    assert not res


@mock.patch("inductiva.users.get_quotas",
            new=fake_get_quotas_no_total_num_machines)
def test__can_start_resource_cant_start_no_num_machines():
    """Test the _can_start_resource function when the quotas dont allow due to
    not enough num machines."""

    resource = mock.MagicMock()
    resource.machine_type = "c2-standard-4"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.04724)
    resource.consumed_vcpus = 4

    res = _can_start_resource(resource)

    assert not res


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_small)
def test__can_start_resource_cant_start_high_cost():
    """Test the _can_start_resource function when the resource cost is too
    high."""

    resource = mock.MagicMock()
    resource.machine_type = "c2-standard-1"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=10.4724)
    resource.consumed_vcpus = 1

    res = _can_start_resource(resource)

    assert not res


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_small)
def test__can_start_resource_can_start_low_quotas():
    """Test the _can_start_resource function when the quotas are small but can
    start machine."""

    resource = mock.MagicMock()
    resource.machine_type = "c2-standard-1"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.04724)
    resource.consumed_vcpus = 1

    res = _can_start_resource(resource)

    assert res


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_normal)
def test__can_start_resource_can_start_too_many_vcpu():
    """Test the _can_start_resource function when the resource asks for
    a valid amount of vcpus."""

    resource = mock.MagicMock()
    resource.machine_type = "c2-standard-4"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.04724)
    resource.consumed_vcpus = 4

    res = _can_start_resource(resource)

    assert res


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_small)
def test__can_start_resource_can_start_high_cost():
    """Test the _can_start_resource function when the resource cost is valid."""

    resource = mock.MagicMock()
    resource.machine_type = "c2-standard-1"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.4724)
    resource.consumed_vcpus = 1

    res = _can_start_resource(resource)

    assert res
