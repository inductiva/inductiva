"""Tests for the Simulator class."""
from typing import Optional

import pytest
from unittest import mock
from pytest import mark

from inductiva import types, simulators, resources
import inductiva


class TesterSimulator(simulators.Simulator):
    """Dummy simulator for testing purposes."""
    pass

def test_override_api_method_prefix():
    simulator = simulators.OpenFOAM()
    assert simulator.api_method_name == \
        "fvm.openfoam_foundation.run_simulation"
    simulator.override_api_method_prefix("windtunnel")
    assert simulator.api_method_name == \
        "windtunnel.openfoam_foundation.run_simulation"


def test_simulator_run__non_mpi_enabled__with_mpi_cluster():
    """Check non-mpi simulator raises error with MPICluster.

    Goal: Verify that simulators without the mpi_enabled decorator raise an
    error stating that MPICluster is not available for this simulator."""

    inductiva.api_key = "dummy"
    cluster = resources.MPICluster(machine_type="c2-standard-16",
                                   num_machines=2, register=False)

    simulator = TesterSimulator()

    with pytest.raises(ValueError) as excinfo:
        simulator.validate_computational_resources(cluster)

    expected_error = "The computational resource is invalid."

    assert expected_error in str(excinfo.value)


def test_simulator_run__mpi_enabled__with_mpi_cluster():
    """Check mpi simulator validates MPICluster.

    Goal: Verify that mpi_enabled simulators correctly validate the MPICluster
    on a simulation.run call."""

    inductiva.api_key = "dummy"
    cluster = resources.MPICluster(machine_type="c2-standard-16",
                       num_machines=2, register=False)

    MPIEnabledSim = simulators.simulator.mpi_enabled(TesterSimulator)
    simulator = MPIEnabledSim()

    try:
        simulator.validate_computational_resources(cluster)
    except:
        assert False, "'validate_computational_resources' raised an exception."



def test_simulator_run__non_mpi_enabled__with_machine_group():
    """Check non-mpi simulator runs correctly with a standard machine group.
    
    Goal: Verify that simulators without the mpi_enabled decorator run normally
    with a standard machine group."""

    inductiva.api_key = "dummy"
    machine_group = resources.MachineGroup(machine_type="c2-standard-16",
                                           num_machines=2, register=False)

    simulator = TesterSimulator()

    try:
        simulator.validate_computational_resources(machine_group)
    except:
        assert False, "'validate_computational_resources' raised an exception."


def test_mpi_enabled__dummy_simulator():
    """Check mpi_enabled decorator works correctly with a dummy Simulator.
    
    Goal: Verify that adding the mpi_enabled decorator to a dummy simulator
    adds a new resource (MPICluster) to the _standard_resources tuple.
    """

    mpi_enabled_sim = simulators.simulator.mpi_enabled(TesterSimulator)

    assert resources.MPICluster in mpi_enabled_sim.get_supported_resources()


@mark.parametrize("simulator", [
    simulators.GROMACS,
    simulators.SplishSplash,
    simulators.FEniCSx,
    simulators.FDS,
    simulators.DualSPHysics,
    simulators.SIMSOPT
])
def test_valid_resources__non_mpi_simulators(simulator):
    """Validate  decorator  in non-MPI simulators.
    
    Goal: Verify that the non MPI-compatible simulators are not decorated with
    the mpi_enabled function and that the _standard_resources only contains
    the standard machines."""
    
    assert not (resources.MPICluster in simulator.get_supported_resources())


@mark.parametrize("simulator", [
    simulators.OpenFOAM,
    simulators.REEF3D,
    simulators.SWASH,
    simulators.XBeach
])
def test_valid_resources__mpi_simulators(simulator):
    """Validate the available machines for MPI simulators.
    
    Goal: Verify that the MPI-compatible simulators are decorated with
    the mpi_enabled function and that the _standard_resources is updated
    correctly."""

    assert resources.MPICluster in simulator.get_supported_resources()
