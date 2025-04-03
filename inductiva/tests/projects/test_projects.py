"""Tests for the projects class"""
import inductiva


def test_create_new_project():
    """Tests if can create new project."""

    my_proj = inductiva.projects.Project("benchmarks-hackathon")

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
        print(proj)
        assert isinstance(proj, inductiva.projects.Project)
