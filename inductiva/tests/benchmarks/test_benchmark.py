"""Tests for the Benchmark class."""
from unittest import mock
import pytest
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
    mocked_benchmark.name = "test_benchmark"
    mocked_benchmark.verbose = False
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
    m8.start.assert_called_once_with(wait_for_quotas=wait_for_quotas,
                                     verbose=False)

    simulator_run_calls = [
        mock.call(input_dir="dir",
                  on=m4,
                  a=1,
                  b=1,
                  project="test_benchmark",
                  resubmit_on_preemption=True,
                  verbose=False),
        mock.call(input_dir="dir",
                  on=m4,
                  a=1,
                  b=4,
                  project="test_benchmark",
                  resubmit_on_preemption=True,
                  verbose=False),
        mock.call(input_dir="dir",
                  on=m8,
                  a=1,
                  b=8,
                  project="test_benchmark",
                  resubmit_on_preemption=True,
                  verbose=False),
    ] * num_repeats
    simulator.run.assert_has_calls(calls=simulator_run_calls, any_order=True)
    assert len(simulator.run.call_args_list) == num_repeats * num_runs


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
