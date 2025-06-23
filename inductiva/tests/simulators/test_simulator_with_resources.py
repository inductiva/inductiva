"""Tests for the Simulator class."""
from argparse import Namespace
from unittest import mock
import inspect
from pathlib import Path
import sys
import uuid

from pytest import mark
import pytest

from inductiva import simulators, resources
import inductiva
from inductiva.client import models


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
        self.simulator = "tester"


def new_machine_init(self, machine_type):
    self.machine_type = machine_type


def test_simulator__wrong_version__raises_error():
    inductiva.set_api_key("dummy")
    with pytest.raises(ValueError) as excinfo:
        inductiva.simulators.CaNS(version="999")
    assert "not available" in str(excinfo.value)


def test_append_simulator_image_suffix_based_on_resource__dev():
    inductiva.set_api_key("dummy")
    #has both cpu and gpu versions
    gmx = inductiva.simulators.GROMACS(use_dev=True)

    mg_gpu = mock.Mock()
    mg_gpu.has_gpu.return_value = True

    mg_no_gpu = mock.Mock()
    mg_no_gpu.has_gpu.return_value = False

    # pylint: disable=protected-access
    sim_image_gpu = gmx._append_simulator_image_suffix_based_on_resource(
        gmx._image_uri, mg_gpu)
    sim_image_no_gpu = gmx._append_simulator_image_suffix_based_on_resource(
        gmx._image_uri, mg_no_gpu)

    assert sim_image_gpu.endswith("_gpu_dev")
    assert sim_image_no_gpu.endswith("_dev")


def test_append_simulator_image_suffix_based_on_resource__not_dev():
    inductiva.set_api_key("dummy")
    #has both cpu and gpu versions
    gmx = inductiva.simulators.GROMACS(use_dev=False)

    mg_gpu = mock.Mock()
    mg_gpu.has_gpu.return_value = True

    # pylint: disable=protected-access
    sim_image_gpu = gmx._append_simulator_image_suffix_based_on_resource(
        gmx._image_uri, mg_gpu)

    assert sim_image_gpu.endswith("_gpu")


@mock.patch("inductiva.resources.MPICluster")
def test_validate_computational_resources__unsupported_resource__raise_error(
        mpi_cluster_mock, list_available_fixture):  # pylint: disable=unused-argument
    """Check non-mpi simulator raises error with MPICluster.

    Goal: Verify that simulators without the mpi_enabled decorator raise an
    error stating that MPICluster is not available for this simulator."""
    inductiva.set_api_key("dummy")
    mock_instance = mpi_cluster_mock()

    simulator = TesterSimulator()

    with pytest.raises(ValueError) as excinfo:
        # pylint: disable=protected-access
        simulator._validate_computational_resources(mock_instance)

    assert "The computational resource is invalid" in str(excinfo.value)


def test_mpi_enabled__dummy_simulator():
    """Check mpi_enabled decorator works correctly with a dummy Simulator.

    Goal: Verify that adding the mpi_enabled decorator to a dummy simulator
    adds a new resource (MPICluster) to the _standard_resources tuple.
    """
    inductiva.set_api_key("dummy")

    mpi_enabled_sim = simulators.simulator.mpi_enabled(TesterSimulator)

    # pylint: disable=protected-access
    assert resources.MPICluster in mpi_enabled_sim.get_supported_resources()


@mark.parametrize("simulator", [
    simulators.GROMACS, simulators.SplishSplash, simulators.FDS,
    simulators.DualSPHysics
])
def test_valid_resources__non_mpi_simulators(simulator):
    """Validate  decorator  in non-MPI simulators.

    Goal: Verify that the non MPI-compatible simulators are not decorated with
    the mpi_enabled function and that the _standard_resources only contains
    the standard machines."""
    inductiva.set_api_key("dummy")

    assert resources.MPICluster not in simulator.get_supported_resources()


def test_validate_computational_resources__none_resource__no_wrapper(
        list_available_fixture):  # pylint: disable=unused-argument
    """Verify that resource cannot be None for a standard simulator."""
    inductiva.set_api_key("dummy")

    with pytest.raises(ValueError) as excinfo:
        simulator = TesterSimulator()
        # pylint: disable=protected-access
        simulator._validate_computational_resources(None)

    assert "The computational resource is invalid" in str(excinfo.value)


def test_validate_computational_resources__none_resource_mpi_wrapped(
        list_available_fixture):  # pylint: disable=unused-argument
    """Verify that resource cannot be None with the mpi_enabled decorator."""
    inductiva.set_api_key("dummy")

    with pytest.raises(ValueError) as excinfo:
        simulator = simulators.simulator.mpi_enabled(TesterSimulator)()
        # pylint: disable=protected-access
        simulator._validate_computational_resources(None)

    assert "The computational resource is invalid" in str(excinfo.value)


def test_validate_computational_resources__valid_machine_group__no_error(
        list_available_fixture):  # pylint: disable=unused-argument
    """Check simulator run correctly with a standard machine group.

    Goal: Verify that simulators with and without the mpi_enabled decorator run
    normally with a standard machine group."""
    inductiva.set_api_key("dummy")

    error_message = "'validate_computational_resources' raised an exception."

    with mock.patch.object(inductiva.resources.MachineGroup, "__init__",
                           new_machine_init):
        machine = inductiva.resources.MachineGroup("c2-standard-16")

        try:
            # pylint: disable=protected-access
            TesterSimulator()._validate_computational_resources(machine)
        except ValueError:
            assert False, error_message

        simulator = simulators.simulator.mpi_enabled(TesterSimulator)

        try:
            # pylint: disable=protected-access
            simulator()._validate_computational_resources(machine)
        except ValueError:
            assert False, error_message


def test_validate_computational_resources__valid_mpi_cluster__no_error(
        list_available_fixture):  # pylint: disable=unused-argument
    """Check mpi-enabled simulator runs correctly with a standard MPICluster.

    Goal: Verify that an mpi simulator correctly validated the MPI Cluster"""
    inductiva.set_api_key("dummy")

    error_message = "'validate_computational_resources' raised an exception."

    with mock.patch.object(inductiva.resources.MPICluster, "__init__",
                           new_machine_init):
        machine = inductiva.resources.MPICluster("c2-standard-16")

        simulator = simulators.simulator.mpi_enabled(TesterSimulator)

        try:
            # pylint: disable=protected-access
            simulator()._validate_computational_resources(machine)
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
    inductiva.set_api_key("dummy")

    assert resources.MPICluster in simulator.get_supported_resources()


@mark.parametrize("resubmit_on_preemption", [None, False, True])
def test_resubmit_on_preemption__is_correctly_handled(resubmit_on_preemption):
    inductiva.set_api_key("dummy")
    # Check that the `resubmit_on_preemption` parameter is present in the
    # `run` method of the simulator and that it is passed correctly to the
    # final api call.

    sim_classes = inspect.getmembers(sys.modules["inductiva.simulators"],
                                     inspect.isclass)

    # Some Mock classes
    class TaskApiMock:
        """Mock class for the TasksApi class."""

        def __init__(self, *_args, **_kwargs):
            pass

        def get_task_status(self, *_args, **_kwargs):
            self._tasks_ahead = 0
            return models.TaskStatus.from_dict({
                "id": "123",
                "status": "pending-input",
                "position_in_queue": {
                    "tasks_ahead": 0
                },
                "is_terminated": True
            })

        def get_task(self, *_args, **_kwargs):
            return Namespace(response=Namespace(
                data=('{"status":"pending-input"}').encode("utf-8")))

    class DefaultDictMock(dict):

        def get(self, *_args):
            return ["1.0.0"]

    resubmit_key = "resubmit_on_preemption"

    mock_mg = mock.Mock()
    mock_mg.id = str(uuid.uuid4())
    mock_mg.has_gpu.return_value = True
    mock_mg.gpu_count.return_value = 1
    mock_mg.available_vcpus = 16

    for sim_name, simcls in sim_classes:

        print(f"Testing simulator: {sim_name}")

        # check that the `resubmit_on_preemption` parameter is present in the
        # `run` method of the simulator
        method_signature = inspect.signature(simcls.run)
        assert resubmit_key in method_signature.parameters

        with mock.patch("inductiva.client") as client_mock, \
            mock.patch("inductiva.simulators.simulator.list_available_images") \
               as list_mock, \
            mock.patch("inductiva.api.methods.submit_request") \
                as submit_mock, \
            mock.patch.object(inductiva.simulators.simulator.Simulator,
                 "_validate_computational_resources",) \
                as validate_resources_mock:

            validate_resources_mock.return_value = None

            client_mock.TasksApi = TaskApiMock
            list_mock.return_value = {"production": DefaultDictMock()}

            submit_mock.return_value = models.TaskSubmittedInfo(
                id="123", status="submitted", is_terminated=False)
            if sim_name == "OpenTelemac":
                sim_obj = simcls(version="1.0.0")
            elif sim_name == "CustomImage":
                sim_obj = simcls(container_image="test")
            else:
                sim_obj = simcls()

            # get positional arguments for the `run` method
            args_spec = inspect.getfullargspec(simcls.run).args
            print(args_spec)
            args = ([],) * (len(args_spec) - 2)  # -2 for self and input_dir

            test_input_dir = Path(__file__).parent / "test_input_dir"
            run_kwargs = {
                "on": mock_mg,
                "remote_assets": [],
            }
            if resubmit_on_preemption is not None:
                run_kwargs[resubmit_key] = resubmit_on_preemption
            if sim_name in ("OpenFOAM", "Delft3D"):
                run_kwargs["commands"] = ["ls"]
            if sim_name in ("CP2K", "OpenSees", "AmrWind"):
                run_kwargs["sim_config_filename"] = "test_config_file"

            # pass remote_assets to coawst to avoid our internal checks
            # that check if the input files are present
            if sim_name == "COAWST":
                run_kwargs["build_coawst_script"] = "hello_world.sh"
                run_kwargs["compile_simulator"] = False

            if sim_name in ("SWAN", "SWASH", "SNLSWAN"):
                args = []
                run_kwargs["sim_config_filename"] = "test_config_file"

            sim_obj.run(test_input_dir, *args, **run_kwargs)

            req_arg = submit_mock.call_args[1]["request"]
            print(req_arg)
            if resubmit_on_preemption is None:
                assert not req_arg.resubmit_on_preemption
            else:
                assert bool(
                    req_arg.resubmit_on_preemption) == resubmit_on_preemption

            submit_mock.assert_called_once()
