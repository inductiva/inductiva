"""Project class"""
import contextvars
import logging
from typing import Optional, Union

import inductiva
from inductiva.client import ApiException
from inductiva.client import models
from inductiva.client.apis.tags import projects_api
from inductiva.client.model.project import Project as ProjectModel

_logger = logging.getLogger(__name__)

_CURRENT_PROJECT = contextvars.ContextVar("current_project", default=None)


def get_current_project():
    """Get the current active project."""
    return _CURRENT_PROJECT.get()


class ProjectInfo:
    """Stores info about the project."""

    def __init__(self, **attrs):
        self.task_by_status = {
            models.TaskStatusCode(attr): int(value)
            for attr, value in attrs.items()
        }


def get_projects():
    """Gets all the user's projects."""
    try:
        _logger.debug("Trying to get remote projects")
        api = projects_api.ProjectsApi(inductiva.api.get_client())
        response = api.get_user_projects()
    except ApiException as ex:
        _logger.error("Failed to get remote projects", exc_info=ex)
        raise ex

    return [Project.from_api_response(resp) for resp in response.body]


class Project:
    """Projects management class.
 
    The Project class allows tasks to be collected under the umbrella
    of a single container object. All tasks submitted in the context
    of an opened project will be associated with that project and can
    be retrieved altogether using the project as a filter. A project
    can be created and managed declaratively or using a context
    manager.
 
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
        self._tasks = None

        model = self._get_model(name)
        if model is None:
            model = self._create_model(name)
        self._update_from_api_response(model)

    @staticmethod
    def from_api_response(model: ProjectModel) -> "Project":
        """Builds a `Project` object from an API response.

        Args:
          model: Info about the project. This should be fetched from the
          backend.

        """
        project = Project.__new__(Project)
        project.append = False
        # pylint: disable=protected-access
        project._token = None
        return project._update_from_api_response(model)
        # pylint: enable=protected-access

    @staticmethod
    def _get_model(name: str) -> ProjectModel:
        """Fetches the project info from the backend.

        Args:
          name: The name of the project

        Returns:
          Project models object with the info about the project.

        """
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
        """Creates a project with the given name on the backend."""
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
        """Updates itself with information fetched from the backend."""
        self._info = ProjectInfo(**model.get("task_status_overview"))
        self.created_at = model.get("created_at")
        self.num_tasks = model.get("num_tasks")
        self.name = model.get("name")
        self.id = model.get("id")

        return self

    def open(self):
        """Opens the project.

        Open the project and make it the active one in the calling
        thread. An opened project will ensure that calls to the
        `get_current_project` function will return this project.
        Consecutive calls to this method are idempotent.

        Raises:
            RuntimeError if another opened project exists, i.e.,
            `get_current_project` returns something other than None.

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

        Calls to the `get_current_project` will return `None` after
        the project is closed.  Consecutive calls to this method are
        idempotent.

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
        """Returns project info.

        It makes no requests to the backend. Instead it uses
        previously stored information. Therefore it can be
        outdated. For a method that requests the most recent data from
        the backend use: `get_info`.

        """
        return self._info

    def get_info(self) -> ProjectInfo:
        """Get project information.

        Get updated information on the project. This method executes a
        call to the backend to retrieve the most recent information
        about this project and stores in internally so that it becomes
        available through the `Project.info` property.  This method is
        suitable when up-to-date information is required.

        """
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

    def get_tasks(self,
                  last_n: int = 5,
                  status: Optional[Union[str, models.TaskStatusCode]] = None):
        """Get the last N submitted tasks to this project.

        Get the last N submitted tasks that belong to this project,
        eventually filtered by status. By default, only the last 5
        submitted tasks are returned, irrespectively of their status.

        Args:
            last_n (int): The number of tasks with repect to the submission
                time to fectch.
            status: Status of the tasks to get. If `None`, tasks with any
                status will be returned.
        """
        return inductiva.tasks.get(last_n=last_n, status=status, project=self)

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.close()
