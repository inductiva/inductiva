"""Tests for the projects class"""
from unittest import mock
import pytest

import inductiva

MOCK_PATH_PROJECTS = "inductiva.projects.project.projects_api.ProjectsApi"
MOCK_PATH_CLIENT = "inductiva.projects.project.inductiva.api.get_client"


def test_open_project_and_close():
    """Tests if the open and close of the project changes the context
    variable properly."""
    with mock.patch(MOCK_PATH_PROJECTS), mock.patch(MOCK_PATH_CLIENT):
        assert inductiva.projects.get_current_project() is None

        with inductiva.projects.Project("test_project", append=True) as project:
            assert inductiva.projects.get_current_project() == project
            assert project.opened

        assert inductiva.projects.get_current_project() is None
        assert not project.opened

        project = inductiva.projects.Project("test_project", append=True)
        project.open()
        assert inductiva.projects.get_current_project() == project
        assert project.opened
        project.close()
        assert inductiva.projects.get_current_project() is None
        assert not project.opened


def test_open_project_without_closing():
    """Tests if opening two projects at the same time fails."""
    with mock.patch(MOCK_PATH_PROJECTS), mock.patch(MOCK_PATH_CLIENT):

        expected_message = "another is active."

        with pytest.raises(Exception) as exc_info:
            project_1 = inductiva.projects.Project(name="p1", append=True)
            project_2 = inductiva.projects.Project(name="p2", append=True)
            project_1.open()
            project_2.open()
        assert expected_message in str(exc_info.value)
        assert inductiva.projects.get_current_project() == project_1

        project_1.close()

        with pytest.raises(Exception) as exc_info:
            with inductiva.projects.Project(name="p1", append=True):
                with inductiva.projects.Project(name="p2", append=True):
                    pass
        assert expected_message in str(exc_info.value)
        assert inductiva.projects.get_current_project() is None


def test_open_project_with_append_raises_exception():
    """Tests if opening a project with append=False raises an exception."""
    with mock.patch(MOCK_PATH_PROJECTS), mock.patch(MOCK_PATH_CLIENT):
        project = inductiva.projects.Project(name="project", append=False)

        with pytest.raises(RuntimeError):
            project.open()

        with pytest.raises(RuntimeError):
            with inductiva.projects.Project(name="project", append=False):
                pass
