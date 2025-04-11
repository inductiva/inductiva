"""Tests for the projects class"""
import inductiva


def test_get_default_project():
    """Tests if can create new project."""

    default_proj = inductiva.projects.Project("default")
    assert isinstance(default_proj, inductiva.projects.Project)
    assert default_proj.name == "default"
    # Make sure the project info can be printed
    # without raising an error
    print(str(default_proj))



def test_get_all_projects():
    """Tests if can get all projects."""

    my_projs = inductiva.projects.get_projects()

    for proj in my_projs:
        # ensure the objects are an instance of the Project class
        assert isinstance(proj, inductiva.projects.Project)
        # ensure that printing the project does not raise an error
        print(proj)


def test_add_task_to_project():
    """Tests if can add task to project."""

    task = inductiva.tasks.get_tasks(project="default", last_n=10)[0]
    previous_project_name = task.info.project
    previous_project = inductiva.projects.Project(previous_project_name)

    new_proj = inductiva.projects.Project("new-project")
    new_proj.add_task(task)

    previous_project.get_tasks(status=task.get_status())
    # Check if the task was added to the new project
    # assert task.info.project == new_proj.name
    # Check if the previous project no longer contains the task
    # assert task.info.project != previous_project.name

    # # Check if the task was added to the project
    # assert task.info.project == my_proj.name
    # # Check if the previous project no longer contains the task
