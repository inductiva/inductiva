"""Tests for the Simulator class."""
import pytest
from pytest import mark
from unittest import mock

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
    """Verify that simulators the mpi_enabled decorator run
    normally with a standard machine group."""

    try:
        simulator = TesterSimulator()
        simulator.validate_computational_resources(None)
    except ValueError:
        assert False, "'validate_computational_resources' raised an exception."


def test_validate_computational_resources__none_resource_mpi_wrapped(
        list_available_fixture):  # pylint: disable=unused-argument
    """Verify that simulators with the mpi_enabled decorator run
    normally with a standard machine group."""

    try:
        simulator = simulators.simulator.mpi_enabled(TesterSimulator)()
        simulator.validate_computational_resources(None)
    except ValueError:
        assert False, "'validate_computational_resources' raised an exception."


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
