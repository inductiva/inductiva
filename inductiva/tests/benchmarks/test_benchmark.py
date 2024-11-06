"""Tests for the Benchmark class."""
from unittest import mock
import pytest
from pathlib import Path
from inductiva.benchmarks import Benchmark
from inductiva.resources import MachineGroup
from inductiva.simulators import Simulator


@pytest.fixture
def benchmark_instance():
    mocked_benchmark = mock.MagicMock(spec=Benchmark)
    mocked_benchmark.runs = []
    mocked_benchmark.simulator = None
    mocked_benchmark.input_dir = None
    mocked_benchmark.on = None
    mocked_benchmark.kwargs = {}
    mocked_benchmark.open = mock.MagicMock(return_value=None)
    return mocked_benchmark


def test_benchmark_single_set_default(benchmark_instance):
    Benchmark.set_default(self=benchmark_instance,
                          simulator="sim",
                          input_dir="dir",
                          on="mg",
                          a=1,
                          b=2)
    assert benchmark_instance.simulator == "sim"
    assert benchmark_instance.input_dir == "dir"
    assert benchmark_instance.on == "mg"
    assert benchmark_instance.kwargs == {"a": 1, "b": 2}


def test_benchmark_partial_set_default(benchmark_instance):
    Benchmark.set_default(self=benchmark_instance, simulator="sim")
    assert benchmark_instance.simulator == "sim"
    assert benchmark_instance.input_dir is None
    assert benchmark_instance.on is None
    assert benchmark_instance.kwargs == {}


def test_benchmark_multiple_set_default(benchmark_instance):
    Benchmark.set_default(self=benchmark_instance, simulator="sim", a=1)
    Benchmark.set_default(self=benchmark_instance, a=1, b=2)
    Benchmark.set_default(self=benchmark_instance, input_dir="dir", on="mg")
    assert benchmark_instance.simulator == "sim"
    assert benchmark_instance.input_dir == "dir"
    assert benchmark_instance.on == "mg"
    assert benchmark_instance.kwargs == {"a": 1, "b": 2}


def test_benchmark_single_add_run(benchmark_instance):
    Benchmark.add_run(self=benchmark_instance,
                      simulator="sim",
                      input_dir="dir",
                      on="mg",
                      a=1,
                      b=2)
    assert benchmark_instance.runs == [("sim", "dir", "mg", {"a": 1, "b": 2})]


def test_benchmark_config_before_add_run(benchmark_instance):
    Benchmark.set_default(self=benchmark_instance,
                          simulator="sim",
                          input_dir="dir",
                          on="mg",
                          a=1,
                          b=2)
    assert benchmark_instance.runs == []
    Benchmark.add_run(self=benchmark_instance)
    assert benchmark_instance.runs == [("sim", "dir", "mg", {"a": 1, "b": 2})]


def test_benchmark_multiple_add_runs(benchmark_instance):
    Benchmark.set_default(self=benchmark_instance,
                          simulator="sim",
                          input_dir="dir",
                          a=1)
    Benchmark.add_run(self=benchmark_instance, on="m2", b=2)
    Benchmark.add_run(self=benchmark_instance, on="m4", b=4)
    assert benchmark_instance.runs == [("sim", "dir", "m2", {
        "a": 1,
        "b": 2
    }), ("sim", "dir", "m4", {
        "a": 1,
        "b": 4
    })]
    Benchmark.set_default(self=benchmark_instance, on="m3")
    Benchmark.add_run(self=benchmark_instance, b=3)
    assert benchmark_instance.runs == [("sim", "dir", "m2", {
        "a": 1,
        "b": 2
    }), ("sim", "dir", "m4", {
        "a": 1,
        "b": 4
    }), ("sim", "dir", "m3", {
        "a": 1,
        "b": 3
    })]


def test_benchmark_run(benchmark_instance):
    simulator = mock.MagicMock(spec=Simulator)
    simulator.run = mock.MagicMock(return_value=None)
    m4 = mock.MagicMock(spec=MachineGroup)
    m8 = mock.MagicMock(spec=MachineGroup)
    m4.start = mock.MagicMock(return_value=None)
    m8.start = mock.MagicMock(return_value=None)
    Benchmark.set_default(self=benchmark_instance,
                          simulator=simulator,
                          input_dir="dir",
                          a=1)
    assert benchmark_instance.runs == []
    Benchmark.add_run(self=benchmark_instance, on=m4, b=4)
    Benchmark.add_run(self=benchmark_instance, on=m8, b=8)
    assert benchmark_instance.runs == [
        (simulator, "dir", m4, {
            "a": 1,
            "b": 4
        }),
        (simulator, "dir", m8, {
            "a": 1,
            "b": 8
        }),
    ]
    Benchmark.run(self=benchmark_instance)
    assert benchmark_instance.runs == []


def test_benchmark_runs_info(benchmark_instance):
    task1 = mock.MagicMock()
    task1.download_inputs = mock.MagicMock(return_value=Path("input_dir_path1"))
    task1.info = mock.MagicMock()
    task1.info.time_metrics.computation_seconds.value = 100
    task1.info.estimated_computation_cost = 10
    task1.info.task_id = "task1"
    task1.info.simulator = "sim1"
    task1.info.executer.vm_type = "vm1"

    task2 = mock.MagicMock()
    task2.download_inputs = mock.MagicMock(return_value=Path("input_dir_path2"))
    task2.info = mock.MagicMock()
    task2.info.time_metrics.computation_seconds.value = 200
    task2.info.estimated_computation_cost = 20
    task2.info.task_id = "task2"
    task2.info.simulator = "sim2"
    task2.info.executer.vm_type = "vm2"

    benchmark_instance.get_tasks = mock.MagicMock(return_value=[task1, task2])

    with mock.patch("builtins.open",
                    mock.mock_open(read_data='{"param": "value"}')):
        info = Benchmark.runs_info(self=benchmark_instance, summary=False)
        assert info == [
            {
                "task_id": "task1",
                "simulator": "sim1",
                "machine_type": "vm1",
                "computation_time": 100,
                "estimated_computation_cost": 10,
                "param": "value",
            },
            {
                "task_id": "task2",
                "simulator": "sim2",
                "machine_type": "vm2",
                "computation_time": 200,
                "estimated_computation_cost": 20,
                "param": "value",
            },
        ]


def test_benchmark_runs_info_summary(benchmark_instance):
    task1 = mock.MagicMock()
    task1.download_inputs = mock.MagicMock(return_value=Path("input_dir_path1"))
    task1.info = mock.MagicMock()
    task1.info.time_metrics.computation_seconds.value = 100
    task1.info.estimated_computation_cost = 10
    task1.info.task_id = "task1"
    task1.info.simulator = "sim1"
    task1.info.executer.vm_type = "vm1"

    task2 = mock.MagicMock()
    task2.download_inputs = mock.MagicMock(return_value=Path("input_dir_path2"))
    task2.info = mock.MagicMock()
    task2.info.time_metrics.computation_seconds.value = 200
    task2.info.estimated_computation_cost = 20
    task2.info.task_id = "task2"
    task2.info.simulator = "sim1"
    task2.info.executer.vm_type = "vm2"

    benchmark_instance.get_tasks = mock.MagicMock(return_value=[task1, task2])

    with mock.patch("builtins.open",
                    mock.mock_open(read_data='{"param": "value"}')):
        info = Benchmark.runs_info(self=benchmark_instance, summary=True)
        assert info == [
            {
                "task_id": "task1",
                "machine_type": "vm1",
                "computation_time": 100,
                "estimated_computation_cost": 10,
            },
            {
                "task_id": "task2",
                "machine_type": "vm2",
                "computation_time": 200,
                "estimated_computation_cost": 20,
            },
        ]
