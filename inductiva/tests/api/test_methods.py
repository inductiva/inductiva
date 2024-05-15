"""Tests for inductiva/api/methods.py"""
from unittest import mock
import pytest

import inductiva


def test_submit_to_closed_project_fails():
    """Test if submitting to a closed project raises an exception."""
    with mock.patch(
            "inductiva.projects.get_current_project") as mock_get_project:
        project = inductiva.projects.Project(name="test_project")
        mock_get_project.return_value = project

        assert not project.opened
        with pytest.raises(RuntimeError):
            inductiva.api.methods.submit_task(None, "dummy_method", {}, None,
                                              "", None, "", "gcp")
