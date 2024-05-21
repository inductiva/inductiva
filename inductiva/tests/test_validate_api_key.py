"""Test for the validate_api_key method."""

import pytest
import inductiva
from unittest.mock import patch
from inductiva.api.methods import validate_api_key


# Test for API key validation
def test_validate_api_key_without_api_key():
    inductiva.set_api_key(None)
    with pytest.raises(ValueError):
        validate_api_key()


# Test for version check on first invocation
def test_version_check_first_invocation():
    inductiva.set_api_key("dummy_api_key")

    # Reset the version_checked flag
    if hasattr(validate_api_key, "version_checked"):
        delattr(validate_api_key, "version_checked")

    patch_path = "inductiva.api.methods.compare_client_and_backend_versions"
    with patch(patch_path) as mock_version_check:
        validate_api_key()
        mock_version_check.assert_called_once()

        # Call again to ensure version check is not performed again
        validate_api_key()
        mock_version_check.assert_called_once()


# Test for correct API configuration return
def test_api_configuration_return():
    inductiva.set_api_key("dummy_api_key")
    config = validate_api_key()
    assert config.api_key["APIKeyHeader"] == "dummy_api_key"
