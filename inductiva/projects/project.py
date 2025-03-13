"""Project class"""
import contextvars
import datetime
import logging
import time
from typing import List, Optional, Union

import inductiva
from inductiva.client import ApiException
from inductiva.client import models
from inductiva.client.apis.paths.tasks import Tasks
from inductiva.client.apis.tags import projects_api
from inductiva.client.model.project import Project as ProjectModel
from inductiva.utils.format_utils import bytes_formatter, currency_formatter, timedelta_formatter

_logger = logging.getLogger(__name__)

_CURRENT_PROJECT = contextvars.ContextVar("current_project", default=None)


def get_current_project():
    """Get the current active project.
    
    Returns the currently active project in the calling thread.  If no
    project has been defined, the default project (`None`) will be
    returned.

    """
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
        """Initialize the Project instance.
        
        Args:
          name (str): The name of the project.
          append (bool): A flag indicating that new tasks can be appended to
              the project. If set to False (default), task submission will fail,
              even though the project might be opened.
        """
        self.append = append
        self._token = None
        self._list_of_tasks = None

        model = self._get_model(name)
        if model is None:
            model = self._create_model(name)
        self._update_from_api_response(model)

    @classmethod
    def from_api_response(cls, model: ProjectModel) -> "Project":
        """Builds a `Project` object from an API response.

        Args:
          model: Info about the project. This should be fetched from the
          backend.

        """
        project = cls.__new__(Project)
        project.append = False
        project._token = None
        return project._update_from_api_response(model)

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
        self.num_tasks = int(model.get("num_tasks"))
        self._name = model.get("name")
        self.id = model.get("id")

        return self

    def open(self):
        """Open the project.

        Open the project and make it the active one in the calling
        thread. An opened project will ensure that calls to the
        `get_current_project` function will return this project.
        Consecutive calls to this method are idempotent.

        Raises:
            RuntimeError if another opened project exists, i.e.,
            `get_current_project` returns something other than None.

        """
        if not self.append:
            raise RuntimeError(
                "Trying to open a project with `append=False`.\n"
                "A Project can only be opened when instantiated with the"
                " `append=True` option.")

        _logger.debug("Opening project %s", self._name)
        current_project = get_current_project()

        if current_project is self:
            _logger.debug("Project is already opened.")
            return

        if current_project is not None:
            raise RuntimeError(
                "Trying to open a project when another is active.")

        self._token = _CURRENT_PROJECT.set(self)

    def close(self):
        """Close the project.

        Calls to the `get_current_project` will return `None` after
        the project is closed.  Consecutive calls to this method are
        idempotent.

        """
        if self._token is None:
            return

        _CURRENT_PROJECT.reset(self._token)
        self._token = None

    @property
    def name(self) -> str:
        return self._name

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
        about this project and stores it internally so that it becomes
        available through the `Project.info` property.  This method is
        suitable when up-to-date information is required.

        """
        name = self._name
        model = self._get_model(name)
        self._update_from_api_response(model)
        return self._info

    def _get_project_cost(self):
        """Returns the estimated project cost.

        The estimated project cost is the combination of the cost of all tasks.
        """
        list_of_tasks = self.list()
        total_cost = 0.0
        for task in list_of_tasks:
            total_cost += task.info.estimated_computation_cost
        return total_cost

    def __str__(self) -> str:
        project_cost = currency_formatter(self._get_project_cost())

        #Count status
        status_counts = {}
        total_files = 0
        total_size = 0

        list_of_tasks = self.list()
        total_duration = datetime.timedelta()
        running_tasks_warning = ""
        for task in list_of_tasks:
            status_counts[task.get_status()] = status_counts.get(
                task.get_status(), 0) + 1
            total_files += task.info.data_metrics.output_total_files.value
            total_size += task.info.data_metrics.output_size_bytes.value
            if task.info.start_time and task.info.end_time:
                start = datetime.datetime.fromisoformat(task.info.start_time)
                end = datetime.datetime.fromisoformat(task.info.end_time)
                duration_task = end - start
                total_duration += duration_task
            else:
                running_tasks_warning = (
                    "Warning: Some tasks may not be finished "
                    "yet. The values presented may change as a result.")

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

        return f"Project '{self._name}' with "\
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
            f"  {k}: {v}" for k, v in self._info.task_by_status.items())
        return header + summary

    def get_tasks(
        self,
        last_n: int = -1,
        force_update=False,
        status: Optional[Union[str,
                               models.TaskStatusCode]] = None) -> List[Tasks]:
        """Get the the tasks of this project.

        Get the tasks that belong to this project,
        eventually filtered by status. By default, it will return all the tasks,
        irrespectively of their status.
        This method will only do a request to the back end if the list of
        tasks is None (never requested the list of tasks) or if `force_update`
        is passed as True.

        Args:
            last_n (int): The number of tasks with repect to the submission
                time to fectch. If `last_n<=0` we fetch all tasks submitted
                to the project.
            force_update (bool): Forces the request to the back end, even if we
                already have a list of tasks associated with this project.
            status: Status of the tasks to get. If `None`, tasks with any
                status will be returned.
        """
        logging.warning("Method `get_tasks` from projects is deprecated.\n"
                        "Use `list` instead.")
        return self.list(last_n=last_n,
                         force_update=force_update,
                         status=status)

    def list(
        self,
        last_n: int = -1,
        force_update=False,
        status: Optional[Union[str,
                               models.TaskStatusCode]] = None) -> List[Tasks]:
        """Get the the tasks of this project.

        Get the tasks that belong to this project,
        eventually filtered by status. By default, it will return all the tasks,
        irrespectively of their status.
        This method will only do a request to the back end if the list of
        tasks is None (never requested the list of tasks) or if `force_update`
        is passed as True.

        Args:
            last_n (int): The number of tasks with repect to the submission
                time to fectch. If `last_n<=0` we fetch all tasks submitted
                to the project.
            force_update (bool): Forces the request to the back end, even if we
                already have a list of tasks associated with this project.
            status: Status of the tasks to get. If `None`, tasks with any
                status will be returned.
        """
        if self._list_of_tasks is None or force_update is True:
            list_of_tasks = inductiva.tasks.get_tasks(last_n=last_n,
                                                      project=self,
                                                      status=status)

            self._list_of_tasks = list_of_tasks
        return self._list_of_tasks

    def wait(self):
        """ Wait for all the tasks in a project to complete."""
        all_tasks = self.list(force_update=True)
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
        list_of_tasks = self.list()

        for task in list_of_tasks:
            task.download_outputs(output_dir=f"{self.name}/{task.id}")

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.close()

    def __eq__(self, other) -> bool:
        return isinstance(other, Project) and \
            self.name == other.name and \
            self.id == other.id
