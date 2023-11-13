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
        mock_instance.compare_client_and_backend_versions.return_value = mock.MagicMock(  # pylint: disable=line-too-long
            body={
                "is_valid": True,
                "server_version": "1.0.0"
            })
        yield mock_instance


# Test successful version match
def test_version_match(mock_api_client):  # pylint: disable=unused-argument,redefined-outer-name
    # No exception should be raised
    methods.compare_client_and_backend_versions("1.0.0")


# Test version mismatch
def test_version_mismatch(mock_api_client):  # pylint: disable=redefined-outer-name
    # Create a mock response object
    mock_response = mock.MagicMock()
    mock_response.status = 406
    mock_response.reason = "Not Acceptable"
    mock_response.getheaders.return_value = {"version": "1.0.1"}

    # Create a mock API response object
    mock_api_response = mock.MagicMock()
    mock_api_response.response = mock_response

    # Create the ApiException with the mock API response
    mock_api_exception = exceptions.ApiException(status=406,
                                                 reason="Not Acceptable",
                                                 api_response=mock_api_response)

    # Set the side effect of the mock API client
    mock_api_client.compare_client_and_backend_versions.side_effect = mock_api_exception  # pylint: disable=line-too-long

    with pytest.raises(RuntimeError) as exc_info:
        methods.compare_client_and_backend_versions("1.0.0")

    assert "Client version 1.0.0 is not compatible" in str(exc_info.value)
    assert "API version 1.0.1" in str(exc_info.value)


def test_invalid_client_version_format(mock_api_client):  # pylint: disable=redefined-outer-name
    mock_api_client.compare_client_and_backend_versions.side_effect = exceptions.ApiException(  # pylint: disable=line-too-long
        status=400, reason="Bad Request")

    with pytest.raises(RuntimeError) as exc_info:
        methods.compare_client_and_backend_versions("invalid-version")

    assert "Bad Request" in str(exc_info.value)
