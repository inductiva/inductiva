"""Tests for the Simulator class."""
from typing import Optional

import pytest
from pytest import mark

from inductiva import types, simulators, resources


def run_test(input_dir: str, on=Optional[types.ComputationalResources]):
    """Test run method to check mpi_disabled decorator."""

    return input_dir


def test_override_api_method_prefix():
    simulator = simulators.OpenFOAM()
    assert simulator.api_method_name == \
        "fvm.openfoam_foundation.run_simulation"
    simulator.override_api_method_prefix("windtunnel")
    assert simulator.api_method_name == \
        "windtunnel.openfoam_foundation.run_simulation"


def test_mpi_disabled__run_test__with_mpi_cluster(mocker):
    """Check that the mpi_disabled decorator raises an error when the
    simulator is not MPI compatible and the user tries to run it on a
    MPICluster."""

    mock_register = mocker.patch.object(resources.MPICluster,
                                        "_register_machine_group",
                                        return_value=("id-resource",
                                                      "name-resource"))

    cluster = resources.MPICluster(machine_type="c2-standard-16",
                                   num_machines=2)

    mpi_disabled__run_test = simulators.simulator.mpi_disabled(run_test)

    with pytest.raises(ValueError) as excinfo:
        mpi_disabled__run_test(input_dir="test", on=cluster)

    assert str(excinfo.value) == "MPI is not available for this simulator. " \
                                "Please use a different computational resource."


def test_mpi_disabled__run_test__with_machine_group(mocker):
    """Check that the mpi_disabled decorator raises an error when the
    simulator is not MPI compatible and the user tries to run it on a
    MPICluster."""

    mock_register = mocker.patch.object(resources.MachineGroup,
                                        "_register_machine_group",
                                        return_value=("id-resource",
                                                      "name-resource"))

    cluster = resources.MachineGroup(machine_type="c2-standard-16",
                                     num_machines=2)

    mpi_disabled__run_test = simulators.simulator.mpi_disabled(run_test)

    input_dir = mpi_disabled__run_test(input_dir="test", on=cluster)

    assert input_dir == "test"


@mark.parametrize("simulator", [
    simulators.GROMACS(),
    simulators.SplishSplash(),
    simulators.FEniCSx(),
    simulators.FDS(),
    simulators.DualSPHysics(),
    simulators.SIMSOPT()
])
def test_mpi_disabled__non_mpi_simulators__with_mpi_cluster(simulator, mocker):
    """Validate mpi_disable decorator is set in non-MPI simulators.
    
    Goal: Verify that the non MPI-compatible simulators are decorated with
    the mpi_disabled function which blocks from using MPICluster. This serves as
    a second layer to check all the simulators are correctly decorated."""

    assert getattr(simulator.run, "mpi_disabled_simulator", False)


@mark.parametrize("simulator", [
    simulators.OpenFOAM(),
    simulators.REEF3D(),
    simulators.SWASH(),
    simulators.XBeach()
])
def test_mpi_disabled__mpi_simulators__with_mpi_cluster(simulator, mocker):
    """Validate mpi_disable decorator is set in non-MPI simulators.
    
    Goal: Verify that the non MPI-compatible simulators are decorated with
    the mpi_disabled function which blocks from using MPICluster. This serves as
    a second layer to check all the simulators are correctly decorated."""

    assert not getattr(simulator.run, "mpi_disabled_simulator", False)
