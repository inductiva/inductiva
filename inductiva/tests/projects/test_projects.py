"""Tests for the projects class"""
import inductiva
import pytest

DEFAULT_PROJECT_NAME = "default"


def test_get_default_project():
    """Tests if can create new project."""

    default_proj = inductiva.projects.Project(DEFAULT_PROJECT_NAME)
    assert isinstance(default_proj, inductiva.projects.Project)
    assert default_proj.name == DEFAULT_PROJECT_NAME
    # Make sure the project info can be printed
    # without raising an error
    assert str(default_proj)


def test_get_all_projects():
    """Tests if can get all projects."""

    my_projs = inductiva.projects.get_projects()

    for proj in my_projs:
        # ensure the objects are an instance of the Project class
        assert isinstance(proj, inductiva.projects.Project)
        # ensure that printing the project does not raise an error
        assert str(proj)


def test_move_task_from_default_project():
    """Tests if can add task to project."""

    prev_proj_tasks = inductiva.tasks.get_tasks(project=DEFAULT_PROJECT_NAME,
                                                status="success")

    # Skip test if no successful tasks are available
    if not prev_proj_tasks:
        pytest.skip("No tasks found with status 'success' for testing")

    task = prev_proj_tasks[0]
    previous_project_name = task.info.project
    previous_project = inductiva.projects.Project(previous_project_name)

    new_proj = inductiva.projects.Project("new-project")
    new_proj.add_task(task)

    # Ask again for the previous project tasks
    prev_proj_tasks = inductiva.tasks.get_tasks(project=DEFAULT_PROJECT_NAME,
                                                status="success")
    # Check if the task was removed from the previous project
    assert not any(t.id == task.id for t in prev_proj_tasks)

    new_proj_tasks = new_proj.get_tasks(status=task.get_status())
    # Check if the task was added to the new project
    assert any(t.id == task.id for t in new_proj_tasks)

    # Move it back to the previous project
    previous_project.add_task(task)
