"""Tests for machine_groups functions"""
import pytest
from unittest import mock
import inductiva.client.models
import inductiva

BASE_RESPONSE = {
    "name": "dummy_name",
    "machine_type": "c2-standard-4",
    "disk_size_gb": 10,
    "id": "dummy_id",
    "creation_timestamp": "2025-06-05T20:16:54.270Z",
    "spot": False,
    "max_vms": 2,
    "num_vms": 2,
    "min_vms": 1,
    "provider_id": "GCP",
    "zone": "europe-west1-b",
    "started": False,
    "machines": [],
}

RESPONSE_1 = {**BASE_RESPONSE, **{"type": "standard", "is_elastic": False}}
RESPONSE_2 = {**BASE_RESPONSE, **{"type": "mpi", "is_elastic": False}}
RESPONSE_3 = {**BASE_RESPONSE, **{"type": "standard", "is_elastic": True}}

RESPONSES = [RESPONSE_1, RESPONSE_2, RESPONSE_3]

EXPECTED_RESULTS = [
    inductiva.resources.MachineGroup, inductiva.resources.MPICluster,
    inductiva.resources.ElasticMachineGroup
]


@pytest.mark.parametrize("response, expected_result",
                         zip(RESPONSES, EXPECTED_RESULTS))
def test_get_by_name(response, expected_result):
    mock_compute_api_path =\
        "inductiva.client.ComputeApi"

    with mock.patch(mock_compute_api_path) as mock_compute_api:
        mock_response = mock.MagicMock()

        mock_response = inductiva.client.models.VMGroupConfig(**response)

        mock_get_vm_group_by_name = mock.MagicMock(return_value=mock_response)

        (mock_compute_api.return_value.get_vm_group_by_name
        ) = mock_get_vm_group_by_name

        with mock.patch(
                "inductiva.resources.machine_groups"
                ".BaseMachineGroup.set_mpi_config") as mock_set_mpi_config:
            mock_set_mpi_config.return_value = inductiva.commands.MPIConfig(
                "4.1.6", np=1, use_hwthread_cpus=True)

            result = inductiva.resources.machine_groups.get_by_name(
                "dummy_name")

    assert isinstance(result, expected_result)


def test_zone_region_mutual_exclusivity():
    """Test that providing both zone and region raises ValueError."""
    with pytest.raises(ValueError,
                       match="Cannot specify both 'zone' and 'region'"):
        inductiva.resources.MachineGroup(machine_type="c2-standard-4",
                                         zone="europe-west1-b",
                                         region="europe-west1")


def test_region_only_valid():
    """Test that providing only region is valid (requires mocking API)."""
    mock_compute_api_path = "inductiva.client.ComputeApi"

    response = {**BASE_RESPONSE, **{"type": "standard", "is_elastic": False}}

    with mock.patch(mock_compute_api_path) as mock_compute_api:
        mock_response = inductiva.client.models.VMGroupConfig(**response)
        mock_register = mock.MagicMock(return_value=mock_response)
        mock_compute_api.return_value.register_vm_group = mock_register
        mock_compute_api.return_value.start_vm_group = mock.MagicMock()

        # Mock the estimate_cloud_cost method to avoid issues with incomplete response
        with mock.patch.object(inductiva.resources.MachineGroup,
                               'estimate_cloud_cost'):
            mg = inductiva.resources.MachineGroup(machine_type="c2-standard-4",
                                                  region="europe-west1")

            # Verify region was passed to register_vm_group
            call_args = mock_register.call_args
            register_request = call_args.kwargs['register_vm_group_request']
            assert register_request.region == "europe-west1"


def test_zone_only_valid():
    """Test that providing only zone is valid (backwards compatibility)."""
    mock_compute_api_path = "inductiva.client.ComputeApi"

    response = {**BASE_RESPONSE, **{"type": "standard", "is_elastic": False}}

    with mock.patch(mock_compute_api_path) as mock_compute_api:
        mock_response = inductiva.client.models.VMGroupConfig(**response)
        mock_register = mock.MagicMock(return_value=mock_response)
        mock_compute_api.return_value.register_vm_group = mock_register
        mock_compute_api.return_value.start_vm_group = mock.MagicMock()

        # Mock the estimate_cloud_cost method to avoid issues with incomplete response
        with mock.patch.object(inductiva.resources.MachineGroup,
                               'estimate_cloud_cost'):
            mg = inductiva.resources.MachineGroup(machine_type="c2-standard-4",
                                                  zone="europe-west1-b")

            # Verify zone was passed to register_vm_group
            call_args = mock_register.call_args
            register_request = call_args.kwargs['register_vm_group_request']
            assert register_request.zone == "europe-west1-b"


def test_byoc_with_region_raises_error():
    """Test that BYOC mode with region raises clear error."""
    with pytest.raises(ValueError,
                       match="BYOC .* mode requires an explicit zone"):
        inductiva.resources.MachineGroup(machine_type="c2-standard-4",
                                         region="europe-west1",
                                         byoc=True)


def test_byoc_with_zone_valid():
    """Test that BYOC mode with zone still works."""
    mock_compute_api_path = "inductiva.client.ComputeApi"
    mock_byoc_path = "inductiva.resources.machine_groups.MachineGroup._register_byoc_gcp"

    with mock.patch(mock_compute_api_path):
        with mock.patch(mock_byoc_path):
            # Mock the available_vcpus property to avoid issues with BYOC initialization
            with mock.patch.object(inductiva.resources.MachineGroup,
                                   'available_vcpus',
                                   new_callable=mock.PropertyMock,
                                   return_value=4):
                mg = inductiva.resources.MachineGroup(
                    machine_type="c2-standard-4",
                    zone="europe-west1-b",
                    byoc=True,
                    provider="GCP")
                # Should not raise any error
                assert mg.zone == "europe-west1-b"
                assert mg.region is None


def test_estimate_machine_cost_with_region():
    """Test estimate_machine_cost with region parameter."""
    mock_compute_api_path = "inductiva.client.ComputeApi"

    with mock.patch(mock_compute_api_path) as mock_compute_api:
        mock_price_response = {"on_demand_price": 0.5, "preemptible_price": 0.1}
        mock_compute_api.return_value.get_instance_price = mock.MagicMock(
            return_value=mock_price_response)

        cost = inductiva.resources.estimate_machine_cost(
            machine_type="c2-standard-4", region="europe-west1")

        assert cost == 0.5

        # Verify region was passed to API
        call_args = mock_compute_api.return_value.get_instance_price.call_args
        assert call_args.kwargs['region'] == "europe-west1"


def test_estimate_machine_cost_zone_region_mutual_exclusivity():
    """Test estimate_machine_cost rejects both zone and region."""
    with pytest.raises(ValueError,
                       match="Cannot specify both 'zone' and 'region'"):
        inductiva.resources.estimate_machine_cost(machine_type="c2-standard-4",
                                                  zone="europe-west1-b",
                                                  region="europe-west1")
