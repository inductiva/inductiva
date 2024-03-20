"""Unit tests for the Computational Resources classes."""
from unittest import mock
import pytest

import inductiva


def fake_register(self, **kwargs):  # pylint: disable = unused-argument
    self._id = "id-resource"  # pylint: disable = protected-access
    self._name = "name-resource"  # pylint: disable = protected-access
    self.register = False


@mock.patch.object(inductiva.resources.MPICluster,
                   attribute="_register_machine_group",
                   new=fake_register)
def test_machines__mpicluster__register():
    """Check the registering of a MPICluster.
    
    Goal: Verify that the MPICluster is initializating and the registration is
    processed, by mocking it instead of calling the API.
    """
    inductiva.api_key = "dummy"
    cluster = inductiva.resources.MPICluster(machine_type="c2-standard-16",
                                             num_machines=2)

    # Check that the cluster has been initialized correctly
    assert cluster.name == "name-resource"  # pylint: disable = protected-access
    assert cluster.id == "id-resource"  # pylint: disable = protected-access
    assert cluster.num_machines == 2
    assert cluster.machine_type == "c2-standard-16"
    assert cluster.register is False


@mock.patch.object(inductiva.resources.MachineGroup,
                   attribute="_register_machine_group",
                   new=fake_register)
def test_machines__machine_group__register():
    """Check the registering of a MachineGroup.
    
    Goal: Verify that the MachineGroup is initializating correctly based on a
    mock registration.
    """

    inductiva.api_key = "dummy"
    cluster = inductiva.resources.MachineGroup(machine_type="c2-standard-16",
                                               num_machines=2)

    # Check that the cluster has been initialized correctly
    assert cluster.name == "name-resource"  # pylint: disable = protected-access
    assert cluster.id == "id-resource"  # pylint: disable = protected-access
    assert cluster.num_machines == 2
    assert cluster.machine_type == "c2-standard-16"
    assert cluster.register is False


@mock.patch.object(inductiva.resources.MachineGroup,
                   attribute="_register_machine_group",
                   new=fake_register)
def test_machines__machine_group__ice_register():
    """Check the registering of a MachineGroup with the ICE provider.
    
    Goal: Verify that the MachineGroup is initializating correctly based on a
    mock registration and the ICE provider.
    """

    inductiva.api_key = "dummy"
    machine = inductiva.resources.MachineGroup(machine_type="c2-highmem-4",
                                               provider="ICE")

    # Check that the cluster has been initialized correctly
    assert machine.name == "name-resource"  # pylint: disable = protected-access
    assert machine.id == "id-resource"  # pylint: disable = protected-access
    assert machine.machine_type == "c2-highmem-4"
    assert machine.register is False
    assert machine.provider == "ICE"


@mock.patch.object(inductiva.resources.MachineGroup,
                   attribute="_register_machine_group",
                   new=fake_register)
def test_machines__ice_register__invalid_args():
    """Check the registering of a MachineGroup with the ICE provider fails with
    num_machines higher than one or spot True.
    
    Goal: Verify that the MachineGroup is initializating correctly based on a
    mock registration and the ICE provider.
    """

    inductiva.api_key = "dummy"
    with pytest.raises(ValueError) as exception:
        inductiva.resources.MachineGroup(machine_type="c2-highmem-4",
                                         num_machines=2,
                                         provider="ICE")

    assert "only supports launching one machine" in str(exception.value)

    with pytest.raises(ValueError) as exception:
        inductiva.resources.MachineGroup(machine_type="c2-highmem-4",
                                         spot=True,
                                         provider="ICE")

    assert "only supports persistent machine" in str(exception.value)
