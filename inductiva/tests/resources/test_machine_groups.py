"""Tests for machine_groups functions"""
import pytest
from unittest import mock
import inductiva

BASE_RESPONSE = {
    "name": "dummy_name",
    "machine_type": "c2-standard-4",
    "disk_size_gb": 10,
    "id": "dummy_id",
    "creation_timestamp": 1234,
    "spot": False,
    "max_vms": 2,
    "num_vms": 2,
    "min_vms": 1,
    "provider_id": "GCP",
}

RESPONSE_1 = {
    **BASE_RESPONSE,
    **{
        "type": "standard",
        "is_elastic": False,
        "provider": "GCP",
    }
}
RESPONSE_2 = {**BASE_RESPONSE, **{"type": "mpi", "is_elastic": False}}
RESPONSE_3 = {**BASE_RESPONSE, **{"type": "standard", "is_elastic": True}}

RESPONSES = [RESPONSE_1, RESPONSE_2, RESPONSE_3]

EXPECTED_RESULTS = [
    inductiva.resources.MachineGroup, inductiva.resources.MPICluster,
    inductiva.resources.ElasticMachineGroup
]


@pytest.mark.parametrize("response, expected_result",
                         zip(RESPONSES, EXPECTED_RESULTS))
@mock.patch("inductiva.client.apis.tags.compute_api.ComputeApi")
def test_get_by_name(mock_compute_api, response, expected_result):
    mock_response = mock.MagicMock()
    mock_response.body = response
    (mock_compute_api.return_value.get_vm_group_by_name.return_value
    ) = mock_response
    result = inductiva.resources.machine_groups.get_by_name("mock_machine_name")

    assert isinstance(result, expected_result)
