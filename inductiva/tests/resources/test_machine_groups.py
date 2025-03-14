"""Tests for machine_groups functions"""
import json
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
        "inductiva.resources.machine_groups.compute_api.ComputeApi"

    with mock.patch(mock_compute_api_path) as mock_compute_api:
        mock_response = mock.MagicMock()

        mock_response.response.data = json.dumps(response,
                                                 indent=2).encode("utf-8")

        mock_get_vm_group_by_name = mock.MagicMock(return_value=mock_response)
        (mock_compute_api.return_value.get_vm_group_by_name
        ) = mock_get_vm_group_by_name

        result = inductiva.resources.machine_groups.get_by_name("dummy_name")

    assert isinstance(result, expected_result)
