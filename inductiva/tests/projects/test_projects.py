"""Tests for the projects class"""
import inductiva


def test_create_new_project():
    """Tests if can create new project."""

    my_proj = inductiva.projects.Project("benchmarks-hackathon")

    # pylint: disable=line-too-long
    expected_str = """Project \'benchmarks-hackathon\' with 8 tasks (id=8b918ccf-84bd-449b-a671-a1bf2ced35ad).

Tasks status:
success: 8

Total number of output files: 64
Total size of output: 11.75 GB

Project duration: 14 minutes and 52 seconds
Project total simulated time: 40 minutes and 10 seconds

Estimated project cost: 1.38 US$
"""
    print(expected_str)
    print(str(my_proj))
    assert str(my_proj) == expected_str


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

    previous_project.get_tasks(status = task.get_status())
    # Check if the task was added to the new project
    # assert task.info.project == new_proj.name
    # Check if the previous project no longer contains the task
    # assert task.info.project != previous_project.name

    # # Check if the task was added to the project
    # assert task.info.project == my_proj.name
    # # Check if the previous project no longer contains the task
