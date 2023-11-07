# pylint: disable=redefined-outer-name,unused-argument
"""Tests for the compare_client_and_backend_versions method."""

import json
import pytest
from unittest.mock import patch, MagicMock

from inductiva.api.methods import compare_client_and_backend_versions


@pytest.fixture
def mock_api_client():
    with patch("inductiva.api.methods.ApiClient") as mock_client:
        mock_instance = MagicMock()
        mock_client.return_value.__enter__.return_value = mock_instance

        mock_instance.call_api.return_value = MagicMock(
            status=200,
            data=json.dumps({
                "is_valid": True,
                "server_version": "1.0.0"
            }).encode("utf-8"))

        yield mock_instance


# Test successful version match
def test_version_match(mock_api_client):
    # No exception should be raised
    compare_client_and_backend_versions("1.0.0")


# Test version mismatch
def test_version_mismatch(mock_api_client):
    mock_api_client.call_api.return_value.status = 200
    mock_api_client.call_api.return_value.data = json.dumps({
        "is_valid": False,
        "server_version": "1.0.1"
    }).encode("utf-8")

    with pytest.raises(RuntimeError) as exc_info:
        compare_client_and_backend_versions("1.0.0")

        assert str(
            exc_info.value).startswith("Client version 1.0.0 is not compatible")


# Test invalid client version format
def test_invalid_client_version_format(mock_api_client):
    mock_api_client.call_api.return_value.status = 400

    with pytest.raises(RuntimeError) as exc_info:
        compare_client_and_backend_versions("invalid-version")

        assert str(exc_info.value).startswith(
            "HTTP error occurred while getting API version.")
