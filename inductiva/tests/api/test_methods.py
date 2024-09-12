"""Tests for inductiva/api/methods.py"""
import uuid
from unittest import mock
import pytest

import inductiva

MOCK_PATH_PROJECTS = "inductiva.projects.project.projects_api.ProjectsApi"
MOCK_PATH_CLIENT = "inductiva.projects.project.inductiva.api.get_client"


def test_submit_to_closed_project_fails():
    """Test if submitting to a closed project raises an exception."""
    with mock.patch("inductiva.projects.get_current_project"
                   ) as mock_get_project, mock.patch(
                       MOCK_PATH_PROJECTS), mock.patch(MOCK_PATH_CLIENT):
        project = inductiva.projects.Project(name="test_project")
        mock_get_project.return_value = project
        mg = mock.Mock()
        mg.id = uuid.uuid4()

        assert not project.opened
        with pytest.raises(RuntimeError):
            inductiva.api.methods.submit_task(None, "dummy_method", {}, mg, "",
                                              None, "", "gcp")
