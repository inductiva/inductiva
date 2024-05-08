"""Project class"""
import contextvars
import logging

import inductiva
from inductiva.client import ApiException
from inductiva.client import models
from inductiva.client.apis.tags import projects_api
from inductiva.client.model.project import Project as ProjectModel

_logger = logging.getLogger(__name__)

_CURRENT_PROJECT = contextvars.ContextVar("current_project", default=None)


def get_current_project():
    """Gets the current project"""
    return _CURRENT_PROJECT.get()


class ProjectInfo:
    """Stores info about the project"""

    def __init__(self, **attrs):
        self.task_by_status = {
            models.TaskStatusCode(attr): int(value)
            for attr, value in attrs.items()
        }


def get_projects():
    """Fteches all the projects of a given user"""
    try:
        _logger.debug("Trying to get remote projects")
        api = projects_api.ProjectsApi(inductiva.api.get_client())
        response = api.get_user_projects()
    except ApiException as ex:
        _logger.error("Failed to get remote projects", exc_info=ex)
        raise ex

    return [Project.from_api_response(resp) for resp in response.body]


class Project:
    """Projects class.
 
    This class manages the cuurent project being used.
 
    Example usage:
 
    >>> project = inductiva.projects.Project("test_project")
    >>> project.start()
    >>> simulator.run(...)  # Submitted to  `'test_project'`
    >>> project.stop()
 
    This is equivalent to:
 
    >>> with inductiva.projects.Project("test_project"):
    >>>    simulator.run(...)
 
    """

    def __init__(self, name: str, *, append: bool = False):
        """
        Args:
          name: Name of the project
          append: If we allow or not new tasks to be logged to the project.
        """
        self.append = append
        self._token = None

        model = self._get_model(name)
        if model is None:
            model = self._create_model(name)
        self._update_from_api_response(model)

    @staticmethod
    def from_api_response(model: ProjectModel) -> "Project":
        project = Project.__new__(Project)
        project.append = False
        # pylint: disable=protected-access
        project._token = None
        return project._update_from_api_response(model)
        # pylint: enable=protected-access

    @staticmethod
    def _get_model(name) -> ProjectModel:
        try:
            _logger.debug("Trying to get remote project %s", name)
            api = projects_api.ProjectsApi(inductiva.api.get_client())
            response = api.get_project({"name": name})
            return response.body
        except ApiException as ex:
            if ex.status != 404:
                _logger.error("Failed to get remote project %s",
                              name,
                              exc_info=ex)
                raise ex
            _logger.debug("Project %s does not exist", name)
        return None

    @staticmethod
    def _create_model(name):
        try:
            _logger.debug("Creating remote project %s", name)
            api = projects_api.ProjectsApi(inductiva.api.get_client())
            response = api.create_project({"name": name})
            return response.body
        except ApiException as ex:
            _logger.error("Failed to create remote project %s",
                          name,
                          exc_info=ex)
            raise RuntimeError(f"Unable to create project {name}") from ex

    def _update_from_api_response(self, model: ProjectModel) -> "Project":
        self._info = ProjectInfo(**model.get("task_status_overview"))
        self.created_at = model.get("created_at")
        self.num_tasks = model.get("num_tasks")
        self.name = model.get("name")
        self.id = model.get("id")

        return self

    def open(self):
        """Opens the project.

        It does this by setting the context var `_CURRENT_PROJECT`.
        """
        _logger.debug("Opening project %s", self.name)
        current_project = get_current_project()

        if current_project is self:
            _logger.debug("Project is already opened.")
            return

        if current_project is not None:
            raise RuntimeError(
                "Trying to open a project when another is active.")

        self._token = _CURRENT_PROJECT.set(self)

    def close(self):
        """Closes the project.

        It does this by reseting the context var `_CURRENT_PROJECT`.
        """
        if self._token is None:
            return

        _CURRENT_PROJECT.reset(self._token)
        self._token = None

    @property
    def opened(self) -> bool:
        """Checks if the project is open."""
        return self._token is not None

    @property
    def info(self) -> ProjectInfo:
        """Returns project info"""
        return self._info

    def get_info(self) -> ProjectInfo:
        """Returns and updates project info."""
        name = self.name
        model = self._get_model(name)
        self._update_from_api_response(model)
        return self._info

    def __str__(self) -> str:
        return f"Project '{self.name}' with "\
               f"{self.num_tasks} tasks (id={self.id})"

    def desc(self) -> str:
        header = str(self) + "\n"
        summary = "\n".join(
            f"  {k}: {v}" for k, v in self._info.task_by_status.items())
        return header + summary

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.close()
