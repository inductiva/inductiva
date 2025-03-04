"""Tests for machine_types functions"""
from unittest import mock
import inductiva
import json

RESPONSE = [{
    "machine_type": "c2-standard-4",
    "num_cpus": 4,
    "ram_gb": 64,
    "price": 0.1,
    "provider_id": "GCP",
    "threads_per_core": 2,
    "region": "europe-west1",
    "spot": False,
    "zone": "europe-west1-b",
    "num_gpus": None,
    "gpu_name": None,
}, {
    "machine_type": "g2-standard-4",
    "num_cpus": 4,
    "ram_gb": 64,
    "price": 0.1,
    "provider_id": "GCP",
    "threads_per_core": 2,
    "region": "europe-west1",
    "spot": False,
    "zone": "europe-west1-b",
    "num_gpus": 1,
    "gpu_name": "nvidia-l4",
}]


def test_get_available_machine_types():
    mock_compute_api_path =\
        "inductiva.resources.utils.compute_api.ComputeApi"
    with mock.patch(mock_compute_api_path) as mock_compute_api:
        mock_inner_response = mock.MagicMock()
        mock_inner_response.data = mock.MagicMock()
        mock_inner_response.data.decode.return_value = json.dumps(RESPONSE)
        mock_response = mock.MagicMock()
        mock_response.response = mock_inner_response
        mock_list_available_machine_types = mock.MagicMock(
            return_value=mock_response)
        (mock_compute_api.return_value.list_available_machine_types
        ) = mock_list_available_machine_types
        result = inductiva.resources.get_available_machine_types(provider="GCP")
        assert result[0]["machine_type"] == "c2-standard-4"
        assert result[0]["num_cpus"] == 4
        assert result[0]["ram_gb"] == 64
        assert result[0]["price"] == 0.1
        assert result[0]["provider_id"] == "GCP"
        assert result[1]["machine_type"] == "g2-standard-4"
        assert result[1]["num_cpus"] == 4
        assert result[1]["ram_gb"] == 64
        assert result[1]["price"] == 0.1
        assert result[1]["provider_id"] == "GCP"
        assert result[1]["num_gpus"] == 1
        assert result[1]["gpu_name"] == "nvidia-l4"
