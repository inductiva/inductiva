"""Project class"""
import logging
import time
from typing import List, Optional

from inductiva import tasks
from inductiva import api as inductiva_api
from inductiva.client.models import TaskStatusCode
from inductiva.client import ApiException
# from inductiva.client import models
from inductiva.client.apis.tags import projects_api
from inductiva.utils.format_utils import currency_formatter

_logger = logging.getLogger(__name__)


def get_projects():
    """Gets all the user's projects."""
    try:
        _logger.debug("Trying to get remote projects")
        api = projects_api.ProjectsApi(inductiva_api.get_client())
        response = api.get_user_projects()
    except ApiException as ex:
        _logger.error("Failed to get remote projects", exc_info=ex)
        raise ex

    # pylint: disable=protected-access
    return [Project._from_api_response(resp) for resp in response.body]


class Project:
    """Projects management class.

    Groups related tasks together under a single project.

    Example:
        project = inductiva.projects.Project("test_project")

        task_1 = simulator.run(...)
        project.add_task(task_1)
        task_2 = simulator.run(...)
        project.add_task(task_2)
    """

    def __init__(self, name: str):
        """Initialize the Project instance.

        Args:
          name (str): The name of the project.
        """
        # If the project already exists, we will load it from the backend.
        self._proj_data = self._get_project(name)
        # Else, we will create a new project.
        if not self._proj_data:
            self._proj_data = self._create_project(name)

    def _get_project(self, name: str):
        """Fetches the project info from the backend."""
        try:
            api = projects_api.ProjectsApi(inductiva_api.get_client())
            response = api.get_project({"name": name})
            print(type(response.body))
            return response.body
        except ApiException as ex:
            if ex.status != 404:
                _logger.error("Failed to get project %s", name, exc_info=ex)
                raise ex
        return None

    def _create_project(self, name):
        """Creates a project with the given name on the backend."""
        try:
            api = projects_api.ProjectsApi(inductiva_api.get_client())
            response = api.create_project({"name": name})
            return response.body
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
        return self._proj_data.get("name")

    @property
    def created_at(self) -> str:
        return self._proj_data.get("created_at")

    @property
    def num_tasks(self) -> int:
        return self._proj_data.get("num_tasks")

    @property
    def id(self) -> str:
        return self._proj_data.get("id")

    @property
    def task_by_status(self) -> dict:
        return {
            TaskStatusCode(attr): int(value) for attr, value in
            self._proj_data.get("task_status_overview").items()
        }

    @property
    def estimated_computation_cost(self) -> float:
        """Returns the estimated project cost.

        Computed as the sum of the estimated computation cost of each task.
        """
        return self._proj_data.get("estimated_computation_cost")

    def __str__(self) -> str:
        formatted_cost = currency_formatter(self.estimated_computation_cost)

        status_report = "\n".join(
            f"  {k}: {v}" for k, v in self.task_by_status.items())

        return f"Project '{self.name}' with "\
               f"{self.num_tasks} tasks (id={self.id}).\n"\
               "\nTasks status:\n"\
               f"{status_report}"\
               f"\nEstimated project cost: {formatted_cost}\n"


#               f"Project total simulated time: {total_duration}\n"\
#               f"\nTotal number of output files: {total_files}\n"\
#               f"Total size of output: {total_size}\n"\
#               f"\nProject duration: {duration}\n"\

    def add_task(self, task: tasks.Task):
        """Adds a task to the project.

        Args:
            task: The task to add to the project.
        """
        try:
            api = projects_api.ProjectsApi(inductiva_api.get_client())
            api.add_task_to_project({"task_id": task.id, "name": self.name})
        except ApiException as ex:
            _logger.error("Failed to add task %s to project %s",
                          task.id,
                          self.name,
                          exc_info=ex)
            raise ex

    def get_tasks(self,
                  last_n: int = -1,
                  status: Optional[str] = None) -> List[tasks.Task]:
        """Get the the tasks of this project.

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
        """ Wait for all the tasks in a project to complete."""
        all_tasks = self.get_tasks()
        print("Waiting for ALL tasks to finish")
        while (not all(x.is_terminal() for x in all_tasks)):
            finished = sum(x.is_terminal() for x in all_tasks)
            print(f"Finished: {finished} Total: {len(all_tasks)}", end="\r")
            time.sleep(5)
        print("All tasks in the project terminated.")

    def download_outputs(self):
        """ Downloads all the outputs for all the tasks in the project.

        All the files will be stored inside
        `inductiva_output/<project_name>/<task_id's>`.
        """
        all_tasks = self.get_tasks()

        for task in all_tasks:
            task.download_outputs(output_dir=f"{self.name}/{task.id}")

    def __eq__(self, other) -> bool:
        return isinstance(other, Project) and \
            self.name == other.name and \
            self.id == other.id
