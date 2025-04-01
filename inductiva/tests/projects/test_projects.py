"""Tests for the projects class"""
import inductiva


def test_create_new_project():
    """Tests if can create new project."""

    my_proj = inductiva.projects.Project("my project")
    assert my_proj.name == "my project"
