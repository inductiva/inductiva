"""Unit tests for the Computational Resources classes."""
from unittest import mock
import pytest

import inductiva
from inductiva.client import models
from inductiva.resources.machine_groups import BaseMachineGroup


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
    self.auto_resize_disk_max_gb = 100
    self._free_space_threshold_gb = 5  # pylint: disable = protected-access
    self._size_increment_gb = 10  # pylint: disable = protected-access


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
    self.auto_resize_disk_max_gb = None
    self._free_space_threshold_gb = 5  # pylint: disable = protected-access
    self._size_increment_gb = 10  # pylint: disable = protected-access


@mock.patch.object(inductiva.resources.MachineGroup, "__init__",
                   fake_init_disk_resizable)
def test_machines__machine_group__dynamic_disk_resize_config__resizable():
    inductiva.set_api_key("dummy")
    machine = inductiva.resources.MachineGroup(machine_type="c2-highmem-4",
                                               provider="GCP")
    #pylint: disable = protected-access
    config = machine._dynamic_disk_resize_config()

    assert "free_space_threshold_gb" in config
    assert "max_disk_size_gb" in config


@mock.patch.object(inductiva.resources.MachineGroup, "__init__",
                   fake_init_disk_not_resizable)
def test_machines__machine_group__dynamic_disk_resize_config__not_resizable():
    inductiva.set_api_key("dummy")
    machine = inductiva.resources.MachineGroup(machine_type="c2-highmem-4",
                                               provider="GCP")
    #pylint: disable = protected-access
    config = machine._dynamic_disk_resize_config()

    assert config is None


def fake_get_quotas_small():
    """Fake get_quotas function."""
    return {
        "max_vcpus":
            models.Quota(
                id="max_vcpus",
                max_allowed=1.0,
                in_use=0.0,
                label="",
                unit="",
                scope=models.QuotaScope.GLOBAL,
            ),
        "max_instances":
            models.Quota(
                id="max_instances",
                max_allowed=1.0,
                in_use=0.0,
                label="",
                unit="",
                scope=models.QuotaScope.GLOBAL,
            ),
        "max_price_hour":
            models.Quota(
                id="max_price_hour",
                max_allowed=1.0,
                in_use=0.0,
                label="",
                unit="",
                scope=models.QuotaScope.GLOBAL,
            ),
        "max_disk_size":
            models.Quota(
                id="max_disk_size",
                max_allowed=1.0,
                in_use=None,
                label="",
                unit="",
                scope=models.QuotaScope.INSTANCE,
            ),
        "mg_max_idle":
            models.Quota(
                id="mg_max_idle",
                max_allowed=1.0,
                in_use=None,
                label="",
                unit="",
                scope=models.QuotaScope.INSTANCE,
            ),
        "mg_max_ttl":
            models.Quota(
                id="mg_max_ttl",
                max_allowed=1.0,
                in_use=None,
                label="",
                unit="",
                scope=models.QuotaScope.INSTANCE,
            ),
    }


def fake_get_quotas_normal():
    """Fake get_quotas function."""
    return {
        "max_vcpus":
            models.Quota(
                id="max_vcpus",
                max_allowed=100.0,
                in_use=0.0,
                label="",
                unit="",
                scope=models.QuotaScope.GLOBAL,
            ),
        "max_instances":
            models.Quota(
                id="max_instances",
                max_allowed=10.0,
                in_use=0.0,
                label="",
                unit="",
                scope=models.QuotaScope.GLOBAL,
            ),
        "max_price_hour":
            models.Quota(
                id="max_price_hour",
                max_allowed=20.0,
                in_use=0.0,
                label="",
                unit="",
                scope=models.QuotaScope.GLOBAL,
            ),
        "max_disk_size":
            models.Quota(
                id="max_disk_size",
                max_allowed=50.0,
                in_use=None,
                label="",
                unit="",
                scope=models.QuotaScope.INSTANCE,
            ),
        "mg_max_idle":
            models.Quota(
                id="mg_max_idle",
                max_allowed=150.0,
                in_use=None,
                label="",
                unit="",
                scope=models.QuotaScope.INSTANCE,
            ),
        "mg_max_ttl":
            models.Quota(
                id="mg_max_ttl",
                max_allowed=48.0,
                in_use=None,
                label="",
                unit="",
                scope=models.QuotaScope.INSTANCE,
            ),
    }


def fake_get_quotas_no_total_num_machines():
    """Fake get_quotas function."""
    return {
        "max_vcpus":
            models.Quota(
                id="max_vcpus",
                max_allowed=100.0,
                in_use=0.0,
                label="",
                unit="",
                scope=models.QuotaScope.GLOBAL,
            ),
        "max_instances":
            models.Quota(
                id="max_instances",
                max_allowed=0.0,
                in_use=0.0,
                label="",
                unit="",
                scope=models.QuotaScope.GLOBAL,
            ),
        "max_price_hour":
            models.Quota(
                id="max_price_hour",
                max_allowed=20.0,
                in_use=0.0,
                label="",
                unit="",
                scope=models.QuotaScope.GLOBAL,
            ),
        "max_disk_size":
            models.Quota(
                id="max_disk_size",
                max_allowed=50.0,
                in_use=None,
                label="",
                unit="",
                scope=models.QuotaScope.INSTANCE,
            ),
        "mg_max_idle":
            models.Quota(
                id="mg_max_idle",
                max_allowed=150.0,
                in_use=None,
                label="",
                unit="",
                scope=models.QuotaScope.INSTANCE,
            ),
        "mg_max_ttl":
            models.Quota(
                id="mg_max_ttl",
                max_allowed=48.0,
                in_use=None,
                label="",
                unit="",
                scope=models.QuotaScope.INSTANCE,
            ),
    }


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_small)
def test_can_start_resource_cant_start_low_quotas():
    """Test the can_start_resource function when the quotas are small."""
    resource = mock.MagicMock(spec=BaseMachineGroup)
    resource.machine_type = "c2-standard-4"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.04724)
    resource.n_vcpus.total = 4
    resource.num_machines = 1

    res = BaseMachineGroup.can_start_resource(resource)

    assert not res


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_normal)
def test_can_start_resource_cant_start_too_many_vcpu():
    """Test the can_start_resource function when the resource asks for
    too many vcpus."""

    resource = mock.MagicMock(spec=BaseMachineGroup)
    resource.machine_type = "c2-standard-400"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.04724)
    resource.n_vcpus.total = 400
    resource.num_machines = 1

    res = BaseMachineGroup.can_start_resource(resource)

    assert not res


@mock.patch("inductiva.users.get_quotas",
            new=fake_get_quotas_no_total_num_machines)
def test_can_start_resource_cant_start_no_num_machines():
    """Test the can_start_resource function when the quotas dont allow due to
    not enough num machines."""

    resource = mock.MagicMock(spec=BaseMachineGroup)
    resource.machine_type = "c2-standard-4"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.04724)
    resource.n_vcpus.total = 4
    resource.num_machines = 1

    res = BaseMachineGroup.can_start_resource(resource)

    assert not res


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_small)
def test_can_start_resource_cant_start_high_cost():
    """Test the can_start_resource function when the resource cost is too
    high."""

    resource = mock.MagicMock(spec=BaseMachineGroup)
    resource.machine_type = "c2-standard-1"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=10.4724)
    resource.n_vcpus.total = 1
    resource.num_machines = 1

    res = BaseMachineGroup.can_start_resource(resource)

    assert not res


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_small)
def test_can_start_resource_can_start_low_quotas():
    """Test the can_start_resource function when the quotas are small but can
    start machine."""

    resource = mock.MagicMock(spec=BaseMachineGroup)
    resource.machine_type = "c2-standard-1"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.04724)
    resource.n_vcpus.total = 1
    resource.num_machines = 1

    res = BaseMachineGroup.can_start_resource(resource)

    assert res


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_normal)
def test_can_start_resource_can_start_too_many_vcpu():
    """Test the can_start_resource function when the resource asks for
    a valid amount of vcpus."""

    resource = mock.MagicMock(spec=BaseMachineGroup)
    resource.machine_type = "c2-standard-4"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.04724)
    resource.n_vcpus.total = 4
    resource.num_machines = 1

    res = BaseMachineGroup.can_start_resource(resource)

    assert res


@mock.patch("inductiva.users.get_quotas", new=fake_get_quotas_small)
def test_can_start_resource_can_start_high_cost():
    """Test the can_start_resource function when the resource cost is valid."""

    resource = mock.MagicMock(spec=BaseMachineGroup)
    resource.machine_type = "c2-standard-1"
    resource.estimate_cloud_cost = mock.MagicMock(return_value=0.4724)
    resource.n_vcpus.total = 1
    resource.num_machines = 1

    res = BaseMachineGroup.can_start_resource(resource)

    assert res
