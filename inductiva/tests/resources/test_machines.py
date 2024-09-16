"""Unit tests for the Computational Resources classes."""
from unittest import mock
import pytest

import inductiva


def fake_init(self, machine_type: str, num_machines=1, provider="GCP"):
    self._id = "id-resource"  # pylint: disable = protected-access
    self._name = "name-resource"  # pylint: disable = protected-access
    self.register = False
    self.machine_type = machine_type
    self.num_machines = num_machines
    self.provider = provider


def fake_register(self, **kwargs):  # pylint: disable = unused-argument
    self._id = "id-resource"  # pylint: disable = protected-access
    self._name = "name-resource"  # pylint: disable = protected-access
    self.register = False


@mock.patch.object(inductiva.resources.MPICluster, "__init__", fake_init)
def test_machines__mpicluster__register():
    """Check the registering of a MPICluster.
    
    Goal: Verify that the MPICluster is initializating and the registration is
    processed, by mocking it instead of calling the API.
    """
    inductiva.set_api_key("dummy")
    cluster = inductiva.resources.MPICluster(machine_type="c2-standard-16",
                                             num_machines=2)

    # Check that the cluster has been initialized correctly
    assert cluster.name == "name-resource"  # pylint: disable = protected-access
    assert cluster.id == "id-resource"  # pylint: disable = protected-access
    assert cluster.num_machines == 2
    assert cluster.machine_type == "c2-standard-16"
    assert cluster.register is False


@mock.patch.object(inductiva.resources.MachineGroup, "__init__", fake_init)
def test_machines__machine_group__register():
    """Check the registering of a MachineGroup.
    
    Goal: Verify that the MachineGroup is initializating correctly based on a
    mock registration.
    """

    inductiva.set_api_key("dummy")
    cluster = inductiva.resources.MachineGroup(machine_type="c2-standard-16",
                                               num_machines=2)

    # Check that the cluster has been initialized correctly
    assert cluster.name == "name-resource"  # pylint: disable = protected-access
    assert cluster.id == "id-resource"  # pylint: disable = protected-access
    assert cluster.num_machines == 2
    assert cluster.machine_type == "c2-standard-16"
    assert cluster.register is False


@mock.patch.object(inductiva.resources.MachineGroup, "__init__", fake_init)
def test_machines__machine_group__ice_register():
    """Check the registering of a MachineGroup with the ICE provider.
    
    Goal: Verify that the MachineGroup is initializating correctly based on a
    mock registration and the ICE provider.
    """

    inductiva.set_api_key("dummy")
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

    inductiva.set_api_key("dummy")

    with pytest.raises(ValueError) as exception:
        inductiva.resources.MachineGroup(machine_type="c2-highmem-4",
                                         spot=True,
                                         provider="ICE")

    assert "only supports persistent machine" in str(exception.value)


@mock.patch.object(inductiva.resources.MachineGroup,
                   attribute="_register_machine_group",
                   new=fake_register)
def test_machines__machine_group__invalid_threads_per_core():
    """Check the registering of a MachineGroup fails with
    threads_per_core different than 1 or 2.
    
    Goal: Verify that the MachineGroup threads_per_core validation is working
    correctly based on a mock registration.
    """

    inductiva.set_api_key("dummy")

    with pytest.raises(ValueError) as exception:
        inductiva.resources.MachineGroup(machine_type="c2-standard-4",
                                         threads_per_core=4)

    assert "`threads_per_core` must be either 1 or 2." in str(exception.value)


def fake_init_disk_resizable(self,
                             machine_type: str,
                             num_machines=1,
                             provider="GCP"):
    self._id = "id-resource"  # pylint: disable = protected-access
    self._name = "name-resource"  # pylint: disable = protected-access
    self.register = False
    self.machine_type = machine_type
    self.num_machines = num_machines
    self.provider = provider
    self.disk_config = inductiva.resources.DiskConfig(max_size_gb=100,
                                                      is_resizable=True)


def fake_init_disk_not_resizable(self,
                                 machine_type: str,
                                 num_machines=1,
                                 provider="GCP"):
    self._id = "id-resource"  # pylint: disable = protected-access
    self._name = "name-resource"  # pylint: disable = protected-access
    self.register = False
    self.machine_type = machine_type
    self.num_machines = num_machines
    self.provider = provider
    self.disk_config = inductiva.resources.DiskConfig(max_size_gb=100,
                                                      is_resizable=False)


@mock.patch.object(inductiva.resources.MachineGroup, "__init__",
                   fake_init_disk_resizable)
def test_machines__machine_group__disk_config_to_dict__resizable():
    inductiva.set_api_key("dummy")
    machine = inductiva.resources.MachineGroup(machine_type="c2-highmem-4",
                                               provider="GCP")

    config = machine._disk_config_to_dict(
    )  # pylint: disable = protected-access

    assert "free_space_threshold_gb" in config
    assert "size_increment_gb" in config
    assert "max_disk_size_gb" in config


@mock.patch.object(inductiva.resources.MachineGroup, "__init__",
                   fake_init_disk_not_resizable)
def test_machines__machine_group__disk_config_to_dict__not_resizable():
    inductiva.set_api_key("dummy")
    machine = inductiva.resources.MachineGroup(machine_type="c2-highmem-4",
                                               provider="GCP")

    config = machine._disk_config_to_dict(
    )  # pylint: disable = protected-access

    assert config is None
