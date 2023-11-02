"""Tests for version checks"""
import pytest
import sys
from io import StringIO
from inductiva.utils.version_check import compare_versions, VersionCheckException


def test_compare_versions_major_too_old():
    with pytest.raises(VersionCheckException) as excinfo:
        compare_versions("1.0.0", "2.0.0")
    assert "The installed version (1.0.0) is too old." in str(excinfo.value)


def test_compare_versions_minor_update_available():
    sys.stdout = captured_output = StringIO()
    compare_versions("1.0.0", "1.1.0")
    output = captured_output.getvalue()
    assert "Update available 1.0.0 -> 1.1.0" in output


def test_compare_versions_patch_update_available():
    sys.stdout = captured_output = StringIO()
    compare_versions("1.0.0", "1.0.1")
    output = captured_output.getvalue()
    assert "Update available 1.0.0 -> 1.0.1" in output


def test_compare_versions_same_version():
    sys.stdout = captured_output = StringIO()
    compare_versions("1.0.0", "1.0.0")
    output = captured_output.getvalue()
    assert output == ""


def test_compare_versions_minor_lower_major_higher():
    sys.stdout = captured_output = StringIO()
    compare_versions("2.0.0", "1.5.0")
    output = captured_output.getvalue()
    assert output == ""
