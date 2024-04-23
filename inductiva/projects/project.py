"""Project class"""
import contextvars

import inductiva
from inductiva.client import ApiException
from inductiva.client.apis.tags import projects_api

CURRENT_PROJECT = contextvars.ContextVar("current_project", default=None)


def get_current_project():
    """Gets the current project"""
    return CURRENT_PROJECT.get()


class Project:
    """Projects class.

    This class manages the cuurent project being used. Which has
    implications to where tasks are submitted.

    Example usage:

    >>> project = inductiva.projects.Project("test_project")
    >>> project.start()
    >>> simulator.run(...)  # Submitted to  `'test_project'`
    >>> project.stop()

    This is equivalent to:

    >>> with inductiva.projects.Project("test_project"):
    >>>    simulator.run(...)

    """

    AVAILABLE_MODES = ["r", "w"]

    def __init__(self, name: str, exists_ok: bool = False, mode: str = "w"):
        """
        Args:
          name: Name of the project
          exists_ok: Controls if we raise an error when a project with the same
          name exists.
          mode: `'r'` we do not allow new tasks to be logged to that project.
          `'w'` we allow tasks to be logged to the project.
        """
        if mode not in self.AVAILABLE_MODES:
            raise ValueError(f"Mode must be in {self.AVAILABLE_MODES}")
        self.name = name
        self.exists_ok = exists_ok
        self.mode = mode

        self._api = projects_api.ProjectsApi(inductiva.api.get_client())
        self._token = None

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.stop()

    def start(self):
        """Starts the project.

        It does this by setting the context var `CURRENT_PROJECT` and,
        if a project does not exists already it will be created.

        """
        if get_current_project() is not None:
            raise RuntimeError(
                "Trying to start a project when another is running.")

        try:
            projects = self._api.get_user_projects().body
        except ApiException as e:
            print("Something went wrong trying to get to your projects")
            raise e

        if self.name in [p["name"] for p in projects]:
            if not self.exists_ok:
                raise ValueError("A project with this name already exists")
            self._token = CURRENT_PROJECT.set(self)
        else:
            try:
                self._api.create_project({"name": self.name})
            except ApiException as e:
                print("Something went wrong while trying to create the project")
                raise e
            self._token = CURRENT_PROJECT.set(self)

    def stop(self):
        """Closes the project.

        It does this by reseting the context var `CURRENT_PROJECT`

        """
        if self._token is None:
            print("You are trying to close a project that is not started.")
            return
        CURRENT_PROJECT.reset(self._token)
        self._token = None
