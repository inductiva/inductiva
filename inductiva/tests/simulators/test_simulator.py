"""Tests for the Simulator class."""
from typing import Optional

import pytest
from unittest import mock
from pytest import mark

from inductiva import types, simulators, resources
import inductiva


class TestSimulator(simulators.Simulator):
    """Dummy simulator class."""

    def __init__(self):
        super().__init__()

    def run(self,
            input_dir: types.Path,
            on: Optional[types.ComputationalResources] = None):

        simulators.simulator.validate_computational_resources(
            on, self._supported_resources)
        return input_dir, on


def test_override_api_method_prefix():
    simulator = simulators.OpenFOAM()
    assert simulator.api_method_name == \
        "fvm.openfoam_foundation.run_simulation"
    simulator.override_api_method_prefix("windtunnel")
    assert simulator.api_method_name == \
        "windtunnel.openfoam_foundation.run_simulation"


@mock.patch.object(resources.MPICluster,
                   "_register_machine_group",
                   return_value=("id-resource", "name-resource"))
def test_simulator_run__non_mpi_enabled__with_mpi_cluster(mocker):  # pylint: disable=unused-argument
    """Check non-mpi simulator raises error with MPICluster.

    Goal: Verify that simulators without the mpi_enabled decorator raise an
    error stating that MPICluster is not available for this simulator."""

    inductiva.api_key = "dummy"
    cluster = resources.MPICluster(machine_type="c2-standard-16",
                                   num_machines=2)

    simulator = TestSimulator()

    with pytest.raises(ValueError) as excinfo:
        simulator.run(input_dir="test", on=cluster)

    expected_error = (f"The computational resource ({cluster}) is not valid "
                      f"for this simulator. Valid computational resources are: "
                      f"{simulator._supported_resources}.")  # pylint: disable=protected-access

    assert str(excinfo.value) == expected_error


@mock.patch.object(resources.MPICluster,
                   "_register_machine_group",
                   return_value=("id-resource", "name-resource"))
def test_simulator_run__mpi_enabled__with_mpi_cluster(mocker):  # pylint: disable=unused-argument
    """Check mpi simulator validates MPICluster.

    Goal: Verify that mpi_enabled simulators correctly validate the MPICluster
    on a simulation.run call."""

    inductiva.api_key = "dummy"
    cluster = resources.MPICluster(machine_type="c2-standard-16",
                                   num_machines=2)

    mpi_enabled_simulator = simulators.simulator.mpi_enabled(TestSimulator)
    simulator = mpi_enabled_simulator()

    input_dir, resource = simulator.run(input_dir="test", on=cluster)

    assert input_dir == "test"
    assert isinstance(resource, resources.MPICluster)


@mock.patch.object(resources.MachineGroup,
                   "_register_machine_group",
                   return_value=("id-resource", "name-resource"))
def test_simulator_run__non_mpi_enabled__with_machine_group(mocker):  # pylint: disable=unused-argument
    """Check non-mpi simulator runs correctly with a standard machine group.
    
    Goal: Verify that simulators without the mpi_enabled decorator run normally
    with a standard machine group."""

    inductiva.api_key = "dummy"
    cluster = resources.MachineGroup(machine_type="c2-standard-16",
                                     num_machines=2)

    simulator = TestSimulator()

    input_dir, resource = simulator.run(input_dir="test", on=cluster)

    assert input_dir == "test"
    assert isinstance(resource, resources.MachineGroup)


def test_mpi_enabled__dummy_simulator():
    """Check mpi_enabled decorator works correctly with a dummy Simulator.
    
    Goal: Verify that adding the mpi_enabled decorator to a dummy simulator
    adds a new resource (MPICluster) to the _standard_resources tuple.
    """

    mpi_enabled_simulator = simulators.simulator.mpi_enabled(TestSimulator)
    expected_supported_resources = {
        resources.MachineGroup, resources.ElasticMachineGroup,
        resources.MPICluster
    }

    assert mpi_enabled_simulator._supported_resources == expected_supported_resources


@mark.parametrize("simulator", [
    simulators.GROMACS(),
    simulators.SplishSplash(),
    simulators.FEniCSx(),
    simulators.FDS(),
    simulators.DualSPHysics(),
    simulators.SIMSOPT()
])
def test_valid_resources__non_mpi_simulators(simulator):
    """Validate  decorator  in non-MPI simulators.
    
    Goal: Verify that the non MPI-compatible simulators are not decorated with
    the mpi_enabled function and that the _standard_resources only contains
    the standard machines."""

    expected_supported_resources = {
        resources.MachineGroup, resources.ElasticMachineGroup
    }

    assert simulator._supported_resources == expected_supported_resources  # pylint: disable=protected-access


@mark.parametrize("simulator", [
    simulators.OpenFOAM(),
    simulators.REEF3D(),
    simulators.SWASH(),
    simulators.XBeach()
])
def test_valid_resources__mpi_simulators(simulator):
    """Validate the available machines for MPI simulators.
    
    Goal: Verify that the MPI-compatible simulators are decorated with
    the mpi_enabled function and that the _standard_resources is updated
    correctly."""

    expected_supported_resources = {
        resources.MachineGroup, resources.ElasticMachineGroup,
        resources.MPICluster
    }

    assert simulator._supported_resources == expected_supported_resources
