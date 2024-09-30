"""Tests for the Simulator class."""
from argparse import Namespace
from unittest import mock
import inspect
from pathlib import Path
import sys
import os
import uuid

from pytest import mark
import pytest

from inductiva import simulators, resources
import inductiva

inductiva.set_api_key("dummy")


@pytest.fixture(name="list_available_fixture")
def _list_available_fixture():
    # Fixture that returns a dictionary with the available images for
    # the TesterSimulator local class.
    ret = {
        "production": {
            "testersimulator": ["1.0.0"]
        },
        "development": {
            "testersimulator": ["1.0.0"]
        }
    }
    fn = "inductiva.simulators.simulator.list_available_images"
    with mock.patch(fn, return_value=ret) as mocker:
        yield mocker


class TesterSimulator(simulators.Simulator):
    """Dummy simulator for testing purposes."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "tester.run_simulation"


def test_override_api_method_prefix():
    simulator = simulators.OpenFOAM()
    assert simulator.api_method_name == \
        "fvm.openfoam_foundation.run_simulation"
    simulator.override_api_method_prefix("windtunnel")
    assert simulator.api_method_name == \
        "windtunnel.openfoam_foundation.run_simulation"


def new_machine_init(self, machine_type):
    self.machine_type = machine_type


@mock.patch("inductiva.resources.MPICluster")
def test_validate_computational_resources__unsupported_resource__raise_error(
        mpi_cluster_mock, list_available_fixture):  # pylint: disable=unused-argument
    """Check non-mpi simulator raises error with MPICluster.

    Goal: Verify that simulators without the mpi_enabled decorator raise an
    error stating that MPICluster is not available for this simulator."""
    mock_instance = mpi_cluster_mock()

    simulator = TesterSimulator()

    with pytest.raises(ValueError) as excinfo:
        simulator.validate_computational_resources(mock_instance)

    assert "The computational resource is invalid" in str(excinfo.value)


def test_mpi_enabled__dummy_simulator():
    """Check mpi_enabled decorator works correctly with a dummy Simulator.

    Goal: Verify that adding the mpi_enabled decorator to a dummy simulator
    adds a new resource (MPICluster) to the _standard_resources tuple.
    """

    mpi_enabled_sim = simulators.simulator.mpi_enabled(TesterSimulator)

    assert resources.MPICluster in mpi_enabled_sim.get_supported_resources()


@mark.parametrize("simulator", [
    simulators.GROMACS, simulators.SplishSplash, simulators.FEniCSx,
    simulators.FDS, simulators.DualSPHysics, simulators.SIMSOPT
])
def test_valid_resources__non_mpi_simulators(simulator):
    """Validate  decorator  in non-MPI simulators.

    Goal: Verify that the non MPI-compatible simulators are not decorated with
    the mpi_enabled function and that the _standard_resources only contains
    the standard machines."""

    assert resources.MPICluster not in simulator.get_supported_resources()


def test_validate_computational_resources__none_resource__no_wrapper(
        list_available_fixture):  # pylint: disable=unused-argument
    """Verify that resource cannot be None for a standard simulator."""

    with pytest.raises(ValueError) as excinfo:
        simulator = TesterSimulator()
        simulator.validate_computational_resources(None)

    assert "The computational resource is invalid" in str(excinfo.value)


def test_validate_computational_resources__none_resource_mpi_wrapped(
        list_available_fixture):  # pylint: disable=unused-argument
    """Verify that resource cannot be None with the mpi_enabled decorator."""

    with pytest.raises(ValueError) as excinfo:
        simulator = simulators.simulator.mpi_enabled(TesterSimulator)()
        simulator.validate_computational_resources(None)

    assert "The computational resource is invalid" in str(excinfo.value)


def test_validate_computational_resources__valid_machine_group__no_error(
        list_available_fixture):  # pylint: disable=unused-argument
    """Check simulator run correctly with a standard machine group.

    Goal: Verify that simulators with and without the mpi_enabled decorator run
    normally with a standard machine group."""

    error_message = "'validate_computational_resources' raised an exception."

    with mock.patch.object(inductiva.resources.MachineGroup, "__init__",
                           new_machine_init):
        machine = inductiva.resources.MachineGroup("c2-standard-16")

        try:
            TesterSimulator().validate_computational_resources(machine)
        except ValueError:
            assert False, error_message

        simulator = simulators.simulator.mpi_enabled(TesterSimulator)

        try:
            simulator().validate_computational_resources(machine)
        except ValueError:
            assert False, error_message


def test_validate_computational_resources__valid_mpi_cluster__no_error(
        list_available_fixture):  # pylint: disable=unused-argument
    """Check mpi-enabled simulator runs correctly with a standard MPICluster.

    Goal: Verify that an mpi simulator correctly validated the MPI Cluster"""

    error_message = "'validate_computational_resources' raised an exception."

    with mock.patch.object(inductiva.resources.MPICluster, "__init__",
                           new_machine_init):
        machine = inductiva.resources.MPICluster("c2-standard-16")

        simulator = simulators.simulator.mpi_enabled(TesterSimulator)

        try:
            simulator().validate_computational_resources(machine)
        except ValueError:
            assert False, error_message


@mark.parametrize("simulator", [
    simulators.OpenFOAM, simulators.REEF3D, simulators.SWASH, simulators.XBeach
])
def test_valid_resources__mpi_simulators(simulator):
    """Validate the available machines for MPI simulators.

    Goal: Verify that the MPI-compatible simulators are decorated with
    the mpi_enabled function and that the _standard_resources is updated
    correctly."""

    assert resources.MPICluster in simulator.get_supported_resources()


@mark.parametrize("resubmit_on_preemption", [None, False, True])
def test_resubmit_on_preemption__is_correctly_handled(resubmit_on_preemption):
    # Check that the `resubmit_on_preemption` parameter is present in the
    # `run` method of the simulator and that it is passed correctly to the
    # final api call.

    sim_classes = inspect.getmembers(sys.modules["inductiva.simulators"],
                                     inspect.isclass)

    # Some Mock classes
    class TaskApiMock:

        def __init__(self, *_args, **_kwargs):
            pass

        def get_task_position_in_queue(self, *_args, **_kwargs):
            return Namespace(body={"tasks_ahead": 0})

    class DefaultDictMock(dict):

        def get(self, *_args):
            return ["1.0.0"]

    resubmit_key = "resubmit_on_preemption"

    mock_mg = mock.Mock()
    mock_mg.id = uuid.uuid4()

    for sim_name, simcls in sim_classes:
        # these 2 classes are not wrappers around the simulators in the backend
        if sim_name in ("FEniCSx", "SIMSOPT"):
            continue

        print(f"Testing simulator: {sim_name}")

        # check that the `resubmit_on_preemption` parameter is present in the
        # `run` method of the simulator
        method_signature = inspect.signature(simcls.run)
        assert resubmit_key in method_signature.parameters

        with mock.patch.dict(os.environ,
                             {"DISABLE_TASK_METADATA_LOGGING": "true"}), \
            mock.patch("inductiva.tasks.task.tasks_api") as taskapi_mock, \
            mock.patch("inductiva.simulators.simulator.list_available_images") \
               as list_mock, \
            mock.patch("inductiva.api.methods.submit_request") \
                as submit_mock, \
            mock.patch.object(inductiva.simulators.simulator.Simulator,
                 "validate_computational_resources",) \
                as validate_resources_mock:

            validate_resources_mock.return_value = None

            taskapi_mock.TasksApi = TaskApiMock
            list_mock.return_value = {"production": DefaultDictMock()}

            submit_mock.return_value = {"id": "123", "status": None}
            if sim_name == "CustomImage":
                sim_obj = simcls(container_image="test")
            else:
                sim_obj = simcls()

            # get positional arguments for the `run` method
            args_spec = inspect.getfullargspec(simcls.run).args
            print(args_spec)
            args = ([],) * (len(args_spec) - 2)  # -2 for self and input_dir

            test_input_dir = Path(__file__).parent / "test_input_dir"
            if resubmit_on_preemption is None:
                # test that the default value of
                # `resubmit_on_preemption` is False
                sim_obj.run(test_input_dir, *args, on=mock_mg)
                req_arg = submit_mock.call_args[1]["request"]
                assert not req_arg[resubmit_key]
            else:
                # test that the value of `resubmit_on_preemption` is passed
                # correctly to the final api call
                sim_obj.run(test_input_dir,
                            *args,
                            on=mock_mg,
                            resubmit_on_preemption=resubmit_on_preemption)
                req_arg = submit_mock.call_args[1]["request"]
                assert bool(req_arg[resubmit_key]) == \
                    resubmit_on_preemption

            submit_mock.assert_called_once()
