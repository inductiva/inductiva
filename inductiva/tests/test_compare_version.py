"""Tests for the compare_client_and_backend_versions method."""

import pytest
from unittest import mock

import inductiva
from inductiva.client import exceptions


@pytest.fixture(name="mock_api_client")
def api_client_fixture():
    """Mock the API client."""
    with mock.patch("inductiva.client.VersionApi") as mock_api:
        mock_instance = mock.MagicMock()
        mock_api.return_value = mock_instance
        mock_instance.compare_client_and_backend_versions.return_value = mock.MagicMock(  # pylint: disable=line-too-long
            body={
                "is_valid": True,
                "server_version": "1.0.0"
            })
        yield mock_instance


# Test successful version match
def test_version_match(mock_api_client):  # pylint: disable=unused-argument
    """Test that the method returns successfully when the versions match."""
    # No exception should be raised
    inductiva.compare_client_and_backend_versions("1.0.0")


# Test version mismatch
def test_version_mismatch(mock_api_client):
    """Method should raise an exception when the versions do not match."""
    # Create a mock response object
    mock_response = mock.MagicMock()
    mock_response.status = 406
    mock_response.reason = "Not Acceptable"

    # Create a mock API response object
    mock_api_response = mock.MagicMock()
    mock_api_response.response = mock_response

    # Create the ApiException with the mock API response
    mock_api_exception = exceptions.ApiException(status=406,
                                                 reason="Not Acceptable",
                                                 http_resp=mock_api_response)

    # Set the side effect of the mock API client
    mock_api_client.compare_client_and_backend_versions.side_effect = mock_api_exception  # pylint: disable=line-too-long

    with pytest.raises(inductiva.VersionError) as exc_info:
        inductiva.compare_client_and_backend_versions("1.0.0")

    assert "inductiva package (version 1.0.0) is outdated" in str(
        exc_info.value)


def test_invalid_client_version_format(mock_api_client):
    """Method should raise an exception when the client version is invalid."""
    mock_api_client.compare_client_and_backend_versions.side_effect = exceptions.ApiException(  # pylint: disable=line-too-long
        status=400, reason="Bad Request")

    with pytest.raises(RuntimeError) as exc_info:
        inductiva.compare_client_and_backend_versions("invalid-version")

    assert "Bad Request" in str(exc_info.value)
