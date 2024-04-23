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

        with inductiva.projects.Project("test_project") as project:
            assert inductiva.projects.get_current_project() == project

        assert inductiva.projects.get_current_project() is None

        project = inductiva.projects.Project("test_project")
        project.start()
        assert inductiva.projects.get_current_project() == project
        project.stop()
        assert inductiva.projects.get_current_project() is None


def test_open_existing_project__exists_ok__false():
    """Tests if opening a project with exists_ok=False raises an
    error"""
    with mock.patch(MOCK_PATH_PROJECTS) as mock_projects, mock.patch(
            MOCK_PATH_CLIENT):
        assert inductiva.projects.get_current_project() is None

        project_name = "test_project"
        mock_response = [{"name": project_name}]
        mock_projects_api = mock.MagicMock()
        mock_projects_api.get_user_projects.return_value.body = mock_response

        mock_projects.return_value = mock_projects_api

        expected_message = "already exists"
        with pytest.raises(ValueError) as exc_info:
            project = inductiva.projects.Project(name=project_name,
                                                 exists_ok=False)
            project.start()
        assert expected_message in str(exc_info.value)

        with pytest.raises(ValueError) as exc_info:
            with inductiva.projects.Project(name=project_name, exists_ok=False):
                pass
        assert expected_message in str(exc_info.value)

        assert inductiva.projects.get_current_project() is None


def test_open_non_existant_project_works():
    """Tests if, even with exists_ok=False, we can open a new project
    with a different name."""
    with mock.patch(MOCK_PATH_PROJECTS) as mock_projects, mock.patch(
            MOCK_PATH_CLIENT):
        assert inductiva.projects.get_current_project() is None
        project_name = "test_project"
        mock_response = [{"name": project_name}]
        mock_projects_api = mock.MagicMock()
        mock_projects_api.get_user_projects.return_value.body = mock_response
        mock_projects.return_value = mock_projects_api

        with inductiva.projects.Project(name=project_name + "_2") as project:
            assert inductiva.projects.get_current_project() == project

        assert inductiva.projects.get_current_project() is None

        project = inductiva.projects.Project(name=project_name + "_2")
        project.start()
        assert inductiva.projects.get_current_project() == project
        project.stop()
        assert inductiva.projects.get_current_project() is None


def test_open_existing_project__exists_ok__true():
    """Tests if opening a project with exists_ok=True works"""
    with mock.patch(MOCK_PATH_PROJECTS) as mock_projects, mock.patch(
            MOCK_PATH_CLIENT):
        assert inductiva.projects.get_current_project() is None

        project_name = "test_project"
        mock_response = [{"name": project_name}]
        mock_projects_api = mock.MagicMock()
        mock_projects_api.get_user_projects.return_value.body = mock_response

        mock_projects.return_value = mock_projects_api

        project = inductiva.projects.Project(name=project_name)
        project.start()
        assert inductiva.projects.get_current_project() == project
        project.stop()

        with inductiva.projects.Project(name=project_name) as project:
            assert inductiva.projects.get_current_project() == project

        assert inductiva.projects.get_current_project() is None


def test_open_project_without_closing():
    """Tests if opening two projects at the same time fails."""
    with mock.patch(MOCK_PATH_PROJECTS), mock.patch(MOCK_PATH_CLIENT):

        expected_message = "Trying to start a project when another is running."

        with pytest.raises(Exception) as exc_info:
            project_1 = inductiva.projects.Project(name="p1")
            project_2 = inductiva.projects.Project(name="p2")
            project_1.start()
            project_2.start()
        assert str(exc_info.value) == expected_message
        assert inductiva.projects.get_current_project() == project_1

        project_1.stop()

        with pytest.raises(Exception) as exc_info:
            with inductiva.projects.Project(name="p1"):
                with inductiva.projects.Project(name="p2"):
                    pass
        assert str(exc_info.value) == expected_message
        assert inductiva.projects.get_current_project() is None
