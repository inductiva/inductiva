"""Project class"""
import datetime
import logging
import time
from typing import List, Optional

from inductiva import tasks
from inductiva import api as inductiva_api
from inductiva.client.models import TaskStatusCode
from inductiva.client import ApiException
# from inductiva.client import models
from inductiva.client.apis.tags import projects_api
from inductiva.utils.format_utils import bytes_formatter, currency_formatter, timedelta_formatter

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

    def _get_project_cost(self):
        """Returns the estimated project cost.

        The estimated project cost is the combination of the cost of all tasks.
        """
        total_cost = 0.0
        for task in self.get_tasks():
            total_cost += task.info.estimated_computation_cost or 0
        return total_cost

    def __str__(self) -> str:
        project_cost = currency_formatter(self._get_project_cost())

        #Count status
        status_counts = {}
        total_files = 0
        total_size = 0
        duration = datetime.timedelta()

        list_of_tasks = self.get_tasks()

        total_duration = datetime.timedelta()
        running_tasks_warning = ""
        for task in list_of_tasks:
            try:
                status_counts[task.get_status()] = status_counts.get(
                    task.get_status(), 0) + 1
                total_files += task.info.data_metrics.output_total_files.value
                total_size += task.info.data_metrics.output_size_bytes.value
                if task.info.start_time and task.info.end_time:
                    start = datetime.datetime.fromisoformat(
                        task.info.start_time)
                    end = datetime.datetime.fromisoformat(task.info.end_time)
                    duration_task = end - start
                    total_duration += duration_task
                else:
                    running_tasks_warning = (
                        "Warning: Some tasks may not be finished "
                        "yet. The values presented may change as a result.")
            # pylint: disable=broad-exception-caught
            # Catch all exceptions to avoid crashing the program
            except Exception as ex:
                _logger.error("Failed to get all the task info for %s",
                              task.id,
                              exc_info=ex)

        if len(list_of_tasks) > 0:
            # get start/end time
            start_project_time = min(
                datetime.datetime.fromisoformat(task.info.start_time)
                for task in list_of_tasks
                if task.info.start_time)
            end_project_time = max(
                datetime.datetime.fromisoformat(task.info.end_time)
                for task in list_of_tasks
                if task.info.end_time)

            duration = end_project_time - start_project_time

        total_size = bytes_formatter(total_size)
        duration = timedelta_formatter(duration)
        total_duration = timedelta_formatter(total_duration)

        status_report = ""
        for status, count in status_counts.items():
            status_report += f"{status}: {count}\n"

        return f"Project '{self.name}' with "\
               f"{self.num_tasks} tasks (id={self.id}).\n"\
               "\nTasks status:\n"\
               f"{status_report}"\
               f"\nTotal number of output files: {total_files}\n"\
               f"Total size of output: {total_size}\n"\
               f"\nProject duration: {duration}\n"\
               f"Project total simulated time: {total_duration}\n"\
               f"\nEstimated project cost: {project_cost}\n"\
               f"{running_tasks_warning}"

    def describe(self) -> str:
        """Generates a string description of the object

        Returns:
          str: A string description of the object. Includes project name,
            total number of tasks and the number of tasks by status.

        """
        header = str(self) + "\n"
        summary = "\n".join(
            f"  {k}: {v}" for k, v in self.task_by_status.items())
        return header + summary

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
