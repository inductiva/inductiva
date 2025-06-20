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
