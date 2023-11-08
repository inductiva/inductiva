# pylint: disable=redefined-outer-name,unused-argument,line-too-long
"""Tests for the compare_client_and_backend_versions method."""

import pytest
from unittest import mock
from inductiva.api import methods
from inductiva.client import exceptions


@pytest.fixture
def mock_api_client():
    with mock.patch("inductiva.api.methods.VersionApi") as mock_api:
        mock_instance = mock.MagicMock()
        mock_api.return_value = mock_instance
        mock_instance.compare_client_and_backend_versions.return_value = mock.MagicMock(
            body={
                "is_valid": True,
                "server_version": "1.0.0"
            })
        yield mock_instance


# Test successful version match
def test_version_match(mock_api_client):
    # No exception should be raised
    methods.compare_client_and_backend_versions("1.0.0")


# Test version mismatch
def test_version_mismatch(mock_api_client):
    mock_api_client.compare_client_and_backend_versions.return_value.body = {
        "is_valid": False,
        "server_version": "1.0.1"
    }

    with pytest.raises(RuntimeError) as exc_info:
        methods.compare_client_and_backend_versions("1.0.0")

        assert str(
            exc_info.value).startswith("Client version 1.0.0 is not compatible")


def test_invalid_client_version_format(mock_api_client):
    mock_api_client.compare_client_and_backend_versions.side_effect = exceptions.ApiException(
        status=400, reason="Bad Request")

    with pytest.raises(RuntimeError) as exc_info:
        methods.compare_client_and_backend_versions("invalid-version")

    assert "Bad Request" in str(exc_info.value)
