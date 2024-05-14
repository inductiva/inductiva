"""Tests for inductiva/api/methods.py"""
from unittest import mock
import pytest

import inductiva


def test_submit_to_project_with_append_false_fails():
    """Tests if submiting a task to a closed project fails"""
    with mock.patch.object(inductiva.projects.Project,
                           "_get_model") as mock_get_model, mock.patch.object(
                               inductiva.projects.Project, "_create_model"):
        mock_get_model.return_value = {
            "task_status_overview": {},
            "created_at": "today",
            "num_tasks": 0,
            "name": "test_project",
            "id": "29034urjeoi49"
        }

        project = inductiva.projects.Project(name="test_project", append=False)

        with project:
            with pytest.raises(PermissionError):
                inductiva.api.methods.submit_task(None, "dummy_method", {},
                                                  None, "", None, "", "gcp")
