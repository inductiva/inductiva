"""Tests for the Benchmark class."""
from unittest import mock
import pytest
from pathlib import Path
from inductiva import simulators, resources
from inductiva.benchmarks import Benchmark


@pytest.fixture(name="benchmark")
def benchmark_fixture():
    mocked_benchmark = mock.MagicMock(spec=Benchmark)
    mocked_benchmark.runs = []
    mocked_benchmark.simulator = None
    mocked_benchmark.input_dir = None
    mocked_benchmark.on = None
    mocked_benchmark.kwargs = {}
    mocked_benchmark.open = mock.MagicMock(return_value=None)
    return mocked_benchmark


def test_benchmark_single_set_default(benchmark):
    Benchmark.set_default(self=benchmark,
                          simulator="sim",
                          input_dir="dir",
                          on="mg",
                          a=1,
                          b=2)
    assert benchmark.simulator == "sim"
    assert benchmark.input_dir == "dir"
    assert benchmark.on == "mg"
    assert benchmark.kwargs == {"a": 1, "b": 2}


def test_benchmark_partial_set_default(benchmark):
    Benchmark.set_default(self=benchmark, simulator="sim")
    assert benchmark.simulator == "sim"
    assert benchmark.input_dir is None
    assert benchmark.on is None
    assert benchmark.kwargs == {}


def test_benchmark_set_default_kwargs(benchmark):
    Benchmark.set_default(self=benchmark, a=1, b=2)
    assert benchmark.kwargs == {"a": 1, "b": 2}
    Benchmark.set_default(self=benchmark, c=3)
    assert benchmark.kwargs == {"a": 1, "b": 2, "c": 3}
    Benchmark.set_default(self=benchmark, c=4)
    assert benchmark.kwargs == {"a": 1, "b": 2, "c": 4}


def test_benchmark_multiple_set_default(benchmark):
    Benchmark.set_default(self=benchmark, simulator="sim", a=1)
    Benchmark.set_default(self=benchmark, a=1, b=2)
    Benchmark.set_default(self=benchmark, input_dir="dir", on="mg")
    assert benchmark.simulator == "sim"
    assert benchmark.input_dir == "dir"
    assert benchmark.on == "mg"
    assert benchmark.kwargs == {"a": 1, "b": 2}


def test_benchmark_single_add_run(benchmark):
    Benchmark.add_run(self=benchmark,
                      simulator="sim",
                      input_dir="dir",
                      on="mg",
                      a=1,
                      b=2)
    assert benchmark.runs == [("sim", "dir", "mg", {"a": 1, "b": 2})]


def test_benchmark_config_before_add_run(benchmark):
    Benchmark.set_default(self=benchmark,
                          simulator="sim",
                          input_dir="dir",
                          on="mg",
                          a=1,
                          b=2)
    assert benchmark.runs == []
    Benchmark.add_run(self=benchmark)
    assert benchmark.runs == [("sim", "dir", "mg", {"a": 1, "b": 2})]


def test_benchmark_multiple_add_runs(benchmark):
    Benchmark.set_default(self=benchmark, simulator="sim", input_dir="dir", a=1)
    Benchmark.add_run(self=benchmark, on="m2", b=2)
    Benchmark.add_run(self=benchmark, on="m4", b=4)
    assert benchmark.runs == [("sim", "dir", "m2", {
        "a": 1,
        "b": 2
    }), ("sim", "dir", "m4", {
        "a": 1,
        "b": 4
    })]
    Benchmark.set_default(self=benchmark, on="m3")
    Benchmark.add_run(self=benchmark, b=3)
    assert benchmark.runs == [("sim", "dir", "m2", {
        "a": 1,
        "b": 2
    }), ("sim", "dir", "m4", {
        "a": 1,
        "b": 4
    }), ("sim", "dir", "m3", {
        "a": 1,
        "b": 3
    })]


@pytest.mark.parametrize("num_repeats, wait_for_quotas", [
    (1, False),
    (5, False),
    (8, True),
    (64, True),
])
def test_benchmark_run(benchmark, num_repeats, wait_for_quotas):
    simulator = mock.MagicMock(spec=simulators.Simulator)
    simulator.run = mock.MagicMock(return_value=None)

    m4 = mock.MagicMock(spec=resources.MachineGroup)
    m4.started = True
    m4.start = mock.MagicMock(return_value=None)

    m8 = mock.MagicMock(spec=resources.MachineGroup)
    m8.started = False
    m8.start = mock.MagicMock(return_value=None)

    Benchmark.set_default(self=benchmark,
                          simulator=simulator,
                          input_dir="dir",
                          a=1)
    assert benchmark.runs == []

    num_runs = 3
    Benchmark.add_run(self=benchmark, on=m4, b=1)
    Benchmark.add_run(self=benchmark, on=m4, b=4)
    Benchmark.add_run(self=benchmark, on=m8, b=8)

    # yapf: disable
    assert benchmark.runs == [(simulator, "dir", m4, {"a": 1, "b": 1}),
                              (simulator, "dir", m4, {"a": 1, "b": 4}),
                              (simulator, "dir", m8, {"a": 1, "b": 8})]
    # yapf: enable

    Benchmark.run(self=benchmark,
                  num_repeats=num_repeats,
                  wait_for_quotas=wait_for_quotas)

    assert benchmark.runs == []

    m4.start.assert_not_called()
    m8.start.assert_called_once_with(wait_for_quotas=wait_for_quotas)

    simulator_run_calls = [
        mock.call(input_dir="dir", on=m4, a=1, b=1),
        mock.call(input_dir="dir", on=m4, a=1, b=4),
        mock.call(input_dir="dir", on=m8, a=1, b=8),
    ] * num_repeats
    simulator.run.assert_has_calls(calls=simulator_run_calls, any_order=True)
    assert len(simulator.run.call_args_list) == num_repeats * num_runs


def test_benchmark_runs_info_select_all(benchmark):
    task1 = mock.MagicMock()
    task1.download_inputs = mock.MagicMock(return_value=Path("input_dir_path1"))
    task1.info = mock.MagicMock()
    task1.info.time_metrics.computation_seconds.value = 100
    task1.info.estimated_computation_cost = 10
    task1.info.task_id = "task1"
    task1.info.simulator = "sim1"
    task1.info.executer.vm_type = "vm1"
    task1.info.to_dict = mock.MagicMock(return_value={"extra_params": {}})

    task2 = mock.MagicMock()
    task2.download_inputs = mock.MagicMock(return_value=Path("input_dir_path2"))
    task2.info = mock.MagicMock()
    task2.info.time_metrics.computation_seconds.value = 200
    task2.info.estimated_computation_cost = 20
    task2.info.task_id = "task2"
    task2.info.simulator = "sim2"
    task2.info.executer.vm_type = "vm2"
    task2.info.to_dict = mock.MagicMock(return_value={"extra_params": {}})

    benchmark.get_tasks = mock.MagicMock(return_value=[task1, task2])

    info = Benchmark.runs_info(self=benchmark, select="all")
    assert info == [
        {
            Benchmark.InfoKey.TASK_ID: "task1",
            Benchmark.InfoKey.SIMULATOR: "sim1",
            Benchmark.InfoKey.MACHINE_TYPE: "vm1",
            Benchmark.InfoKey.TIME: 100,
            Benchmark.InfoKey.COST: 10,
            "param": "value",
        },
        {
            Benchmark.InfoKey.TASK_ID: "task2",
            Benchmark.InfoKey.SIMULATOR: "sim2",
            Benchmark.InfoKey.MACHINE_TYPE: "vm2",
            Benchmark.InfoKey.TIME: 200,
            Benchmark.InfoKey.COST: 20,
            "param": "value",
        },
    ]


def test_benchmark_runs_info_select_distinct(benchmark):
    task1 = mock.MagicMock()
    task1.download_inputs = mock.MagicMock(return_value=Path("input_dir_path1"))
    task1.info = mock.MagicMock()
    task1.info.time_metrics.computation_seconds.value = 100
    task1.info.estimated_computation_cost = 10
    task1.info.task_id = "task1"
    task1.info.simulator = "sim1"
    task1.info.executer.vm_type = "vm1"
    task1.info.to_dict = mock.MagicMock(return_value={"extra_params": {}})

    task2 = mock.MagicMock()
    task2.download_inputs = mock.MagicMock(return_value=Path("input_dir_path2"))
    task2.info = mock.MagicMock()
    task2.info.time_metrics.computation_seconds.value = 200
    task2.info.estimated_computation_cost = 20
    task2.info.task_id = "task2"
    task2.info.simulator = "sim1"
    task2.info.executer.vm_type = "vm2"
    task2.info.to_dict = mock.MagicMock(return_value={"extra_params": {}})

    benchmark.get_tasks = mock.MagicMock(return_value=[task1, task2])

    info = Benchmark.runs_info(self=benchmark, select="distinct")
    assert info == [
        {
            Benchmark.InfoKey.TASK_ID: "task1",
            Benchmark.InfoKey.MACHINE_TYPE: "vm1",
            Benchmark.InfoKey.TIME: 100,
            Benchmark.InfoKey.COST: 10,
        },
        {
            Benchmark.InfoKey.TASK_ID: "task2",
            Benchmark.InfoKey.MACHINE_TYPE: "vm2",
            Benchmark.InfoKey.TIME: 200,
            Benchmark.InfoKey.COST: 20,
        },
    ]


def test_benchmark_terminate(benchmark):
    task1 = mock.MagicMock()
    task1.info = mock.MagicMock()
    task1.info.executer.vm_name = "vm1"

    task2 = mock.MagicMock()
    task2.info = mock.MagicMock()
    task2.info.executer.vm_name = "vm2"

    benchmark.get_tasks = mock.MagicMock(return_value=[task1, task2])

    machine1 = mock.MagicMock()
    machine1.name = "vm1"
    machine1.terminate = mock.MagicMock()

    machine2 = mock.MagicMock()
    machine2.name = "vm2"
    machine2.terminate = mock.MagicMock()

    resources.get = mock.MagicMock(return_value=[machine1, machine2])

    Benchmark.terminate(self=benchmark)

    machine1.terminate.assert_called_once_with(verbose=False)
    machine2.terminate.assert_called_once_with(verbose=False)


def test_benchmark_terminate_no_tasks(benchmark):
    benchmark.get_tasks = mock.MagicMock(return_value=[])

    machine1 = mock.MagicMock()
    machine1.name = "vm1"
    machine1.terminate = mock.MagicMock()

    machine2 = mock.MagicMock()
    machine2.name = "vm2"
    machine2.terminate = mock.MagicMock()

    resources.get = mock.MagicMock(return_value=[machine1, machine2])

    Benchmark.terminate(self=benchmark)

    machine1.terminate.assert_not_called()
    machine2.terminate.assert_not_called()


def test_benchmark_terminate_single_terminated_task(benchmark):
    task = mock.MagicMock()
    task.info = mock.MagicMock()
    task.info.executer.vm_name = "vm1"

    benchmark.get_tasks = mock.MagicMock(return_value=[task])

    machine = mock.MagicMock()
    machine.name = "vm2"
    machine.terminate = mock.MagicMock()

    resources.get = mock.MagicMock(return_value=[machine])

    Benchmark.terminate(self=benchmark)

    machine.terminate.assert_not_called()


def test_benchmark_terminate_duplicate_tasks(benchmark):
    task1 = mock.MagicMock()
    task1.info = mock.MagicMock()
    task1.info.executer.vm_name = "vm1"

    task2 = mock.MagicMock()
    task2.info = mock.MagicMock()
    task2.info.executer.vm_name = "vm1"

    benchmark.get_tasks = mock.MagicMock(return_value=[task1, task2])

    machine = mock.MagicMock()
    machine.name = "vm1"
    machine.terminate = mock.MagicMock()

    resources.get = mock.MagicMock(return_value=[machine])

    Benchmark.terminate(self=benchmark)

    machine.terminate.assert_called_once_with(verbose=False)
