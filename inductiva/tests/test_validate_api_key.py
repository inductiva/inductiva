"""Test for the get_api_config method."""

import pytest
import inductiva
from unittest.mock import patch
from inductiva.api.methods import get_api_config
from inductiva import ApiKeyError, _set_key_and_check_version


# Test for API key validation
def test_set_api_key_config_without_api_key():
    with pytest.raises(ApiKeyError):
        inductiva.set_api_key(None)


# Test for version check on first invocation
def test_version_check_first_invocation():
    inductiva.set_api_key("dummy_api_key")
    # Reset the version_checked flag
    if hasattr(_set_key_and_check_version, "version_checked"):
        delattr(_set_key_and_check_version, "version_checked")

    patch_path = "inductiva.compare_client_and_backend_versions"
    with patch(patch_path) as mock_version_check:
        _set_key_and_check_version()
        mock_version_check.assert_called_once()

        # Call again to ensure version check is not performed again
        _set_key_and_check_version()
        mock_version_check.assert_called_once()


# Test for correct API configuration return
def test_api_configuration_return():
    inductiva.set_api_key("dummy_api_key")
    config = get_api_config()
    assert config.api_key["APIKeyHeader"] == "dummy_api_key"
