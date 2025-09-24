"""
This module provides functionality for managing projects and their associated
tasks within the Inductiva platform. A project serves as a container for
grouping related tasks, enabling better organization and management of
computational workflows.

Classes:
    - Project: Represents a project that groups related tasks together.

Functions:
    - get_projects(): Retrieves all projects associated with the current user.

Key Features:
    - Create and manage projects on the backend.
    - Add tasks to projects and retrieve the most recent tasks or filter by
      status.
    - Monitor the status of all tasks in a project and wait for their
    completion.
    - Download outputs for all tasks in a project.
    - Estimate the total computation cost of a project based on its tasks.

Example Usage:

      .. code-block:: python

        import inductiva

        # Create a new project or load an existing one
        project = inductiva.projects.Project("my_project")

        # Select a simulator
        fvcom = inductiva.simulators.FVCOM()
        # Run and add tasks to the project
        for i in range(10):
            task = fvcom.run(...)
            project.add_task(task)

        # Monitor task completion
        project.wait()
        # Download outputs for all tasks
        project.download_outputs()
        # Print project details
        print(project)
"""

import logging
import time
from typing import List, Optional

from inductiva import tasks
from inductiva import api as inductiva_api
import inductiva.client
from inductiva.client.models import TaskStatusCode, ProjectCreate, ProjectType
from inductiva.client import ApiException

# from inductiva.client import models
from inductiva.utils import files, format_utils

_logger = logging.getLogger(__name__)


def get_projects() -> List["Project"]:
    """Gets all the user's projects."""
    try:
        _logger.debug("Trying to get remote projects")
        api = inductiva.client.ProjectsApi(inductiva_api.get_client())
        response = api.get_user_projects()
    except ApiException as ex:
        _logger.error("Failed to get remote projects", exc_info=ex)
        raise ex

    # pylint: disable=protected-access
    return [Project._from_api_response(resp) for resp in response]


class Project:
    """
    Projects management class.

    Groups related tasks together under a single project.

    Example:

      .. code-block:: python

        project = inductiva.projects.Project("my project")

        task_1 = simulator.run(...)
        project.add_task(task_1)

        task_2 = simulator.run(...)
        project.add_task(task_2)
    """

    def __init__(self, name: str):
        """
        Initialize the Project instance.

        Args:
          name (str): The name of the project.
        """
        self._api = inductiva.client.ProjectsApi(inductiva_api.get_client())
        # If the project already exists, we will load it from the backend.
        self._proj_data = self._get_project(name)
        # Else, we will create a new project.
        if not self._proj_data:
            self._proj_data = self._create_project(name)

    def _get_project_type(self):
        return ProjectType.PROJECT

    def _get_project(self, name: str):
        """Fetches the project info from the backend."""
        try:
            return self._api.get_project(name=name)
        except ApiException as ex:
            if ex.status != 404:
                _logger.error("Failed to get project %s", name, exc_info=ex)
                raise ex
        return None

    def _create_project(self, name):
        """Creates a project with the given name on the backend."""
        try:
            return self._api.create_project(project_create=ProjectCreate(
                name=name, project_type=self._get_project_type()))
        except ApiException as ex:
            _logger.error("Failed to create project %s", name, exc_info=ex)
            raise RuntimeError(f"Unable to create project {name}") from ex

    @classmethod
    def _from_api_response(cls, resp):
        """Creates a Project instance from the API response."""
        project = cls.__new__(cls)
        project._proj_data = resp
        return project

    @property
    def name(self) -> str:
        """Returns the name of the project."""
        return self._proj_data.name

    @property
    def created_at(self) -> str:
        """Returns the creation date and time of the project."""
        return self._proj_data.created_at.isoformat()

    @property
    def num_tasks(self) -> int:
        """Returns the number of tasks in the project."""
        return self._proj_data.num_tasks

    @property
    def id(self) -> str:
        """Returns the unique ID of the project."""
        return self._proj_data.id

    @property
    def task_by_status(self) -> dict:
        """
        Returns a dictionary with the number of tasks by status.
        The keys are the status codes and the values are the number of tasks
        with that status.
        """
        return {
            TaskStatusCode(attr): int(value)
            for attr, value in self._proj_data.task_status_overview.items()
        }

    @property
    def estimated_computation_cost(self) -> float:
        """
        Returns the estimated project cost.

        Computed as the sum of the estimated computation cost of each task.
        """
        return self._proj_data.estimated_computation_cost

    def __str__(self) -> str:
        formatted_cost = format_utils.currency_formatter(
            self.estimated_computation_cost)
        formatted_created_at = format_utils.datetime_formatter_ymd_hm(
            self.created_at)
        status_report = "\n".join(
            f"  {k}: {v}" for k, v in self.task_by_status.items())

        return (f"Project '{self.name}' created at {formatted_created_at}.\n"
                f"\nTotal number of tasks: {self.num_tasks}\n"
                "\nTasks by status:\n"
                f"{status_report}\n"
                f"\nEstimated total computation cost: {formatted_cost}\n")

    def add_task(self, task: tasks.Task):
        """
        Adds a task to the project.

        Args:
            task: The task to add to the project.
        """
        try:
            self._api.add_task_to_project(name=self.name, task_id=task.id)
        except ApiException as ex:
            _logger.error(
                "Failed to add task %s to project %s",
                task.id,
                self.name,
                exc_info=ex,
            )
            raise ex

    def get_tasks(self,
                  last_n: int = -1,
                  status: Optional[str] = None) -> List[tasks.Task]:
        """
        Get the the tasks of this project.

        Optionally, those can be filtered by task status.

        Args:
            last_n (int): The number of tasks with repect to the submission
                time to fetch. If `last_n<=0` we fetch all tasks submitted
                to the project.
            status: Status of the tasks to get. If `None`, tasks with any
                status will be returned.
        """
        return tasks.get_tasks(last_n=last_n, project=self.name, status=status)

    def wait(self):
        """Wait for all the tasks in a project to complete."""
        all_tasks = self.get_tasks()
        print("Waiting for ALL tasks to finish")
        while not all(x.is_terminal() for x in all_tasks):
            finished = sum(x.is_terminal() for x in all_tasks)
            print(f"Finished: {finished} Total: {len(all_tasks)}", end="\r")
            time.sleep(5)
        print("All tasks in the project terminated.")

    def download_outputs(self, output_dir: Optional[str] = None):
        """Downloads all the outputs for all the tasks in the project.

        All task outputs will be organized within the specified `output_dir`.
        If `output_dir` is not provided, outputs will be saved to a default
        location under `inductiva_output/<project_name>/<task_id>/`.
        Otherwise, they will be stored in `<output_dir>/<task_id>/`.

        Args:
            output_dir (str, optional): The base directory where project outputs
                                        will be downloaded.
        """

        for task in self.get_tasks():
            base_path = output_dir or files.resolve_output_path(self.name)
            task.download_outputs(output_dir=f"{base_path}/{task.id}")

    def delete(self):
        """Delete a project on the backend.
        
        This method does not delete the project tasks, only the project itself.
        The tasks will be moved to the "default" project.
        """
        try:
            return self._api.delete_project(project_name=self.name)
        except ApiException as ex:
            _logger.error("Failed to delete project %s", self.name, exc_info=ex)
            raise RuntimeError(f"Unable to delete project {self.name}") from ex

    def __eq__(self, other) -> bool:
        return (isinstance(other, Project) and self.name == other.name and
                self.id == other.id)
