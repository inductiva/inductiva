"""Project class"""
import contextvars
import logging

import inductiva
from inductiva import tasks
from inductiva.client import models
from inductiva.client import ApiException
from inductiva.client.apis.tags import projects_api
from inductiva.client.model.project import Project as ProjectModel


_logger = logging.getLogger(__name__)

_CURRENT_PROJECT = contextvars.ContextVar("current_project", default=None)


def get_current_project():
    """Gets the current project"""
    return _CURRENT_PROJECT.get()


class ProjectInfo:
    def __init__(self, **attrs):
        self.task_by_status = {
            models.TaskStatusCode(attr): int(value)
            for attr, value in attrs.items()
        }


def get_projects():
    try:
        _logger.debug("Trying to get remote projects")
        api = projects_api.ProjectsApi(inductiva.api.get_client())
        response = api.get_user_projects()
    except ApiException as ex:
        _logger.error("Failed to get remote projects", exc_info=ex)
        raise ex
    
    return [Project.from_api_response(resp) for resp in response.body]


class Project:
    def __init__(self, name: str, *, append: bool = False):
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
        project._token = None
        return project._update_from_api_response(model)

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
                               name, exc_info=ex)
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
                          name, exc_info=ex)
            raise RuntimeError(f"Unable to create project {name}") from ex

    def _update_from_api_response(self, model: ProjectModel) -> "Project":
        self._info = ProjectInfo(**model.get('task_status_overview'))
        self.created_at = model.get('created_at')
        self.num_tasks = model.get('num_tasks')
        self.name = model.get('name')
        self.id = model.get('id')

        return self

    def open(self):
        _logger.debug("Opening project %s", self.name)
        current_project = get_current_project()

        if current_project is self:
            _logger.debug("Project is already opened.")
            return

        if current_project is not None:
            raise RuntimeError("Trying to open a project when another is active.")

        self._token = _CURRENT_PROJECT.set(self)


    def close(self):
        if self._token is None:
            return

        _CURRENT_PROJECT.reset(self._token)
        self._token = None

    @property
    def opened(self) -> bool:
        # current_project = get_current_project()
        # if current_project is None:
        #     return False
        # return current_project.id == self.id
        return self._token is not None

    @property
    def info(self) -> ProjectInfo:
        return self._info
    
    def get_info(self) -> ProjectInfo:
        self._info = ...
        return self._info

    # @property
    # def tasks(self):
    #     if self._tasks is None:
    #         return self.get_tasks()
    #     return self._tasks
    
    # def get_tasks(self) -> List[Task]:
    #    tasks
       
    #    self._tasks = tasks
    #    return self._tasks


    def __str__(self) -> str:
        return f"Project '{self.name}' with {self.num_tasks} tasks (id={self.id})"

    def desc(self) -> str:
        header = str(self) + "\n"
        summary = "\n".join(f"  {k}: {v}" for k, v in self._info.task_by_status.items())
        return header + summary
    
    def __enter__(self):
        self.open()
        return self
    
    def __exit__(self, exc_type, exc_value, exc_tb):
        self.close()

_logger.setLevel(logging.DEBUG)
import sys
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)




# p = Project("ssantosa716d17d")#, exists_ok=False)
# print(p)
# print(p.info)
# print(p)
# print(get_current_project())

# print(p.desc())
# print(p.opened)
# p.open()
# print(p.opened)




p = Project("ssantosa716d17d")
print(p)

projects = get_projects()
print(projects)

p.open()
print(p.opened)

for x in projects:
    print(x)
    print(x.opened)



# print(p)
# print(f"{p.opened=}")
# print(f"{get_current_project()=}")
# with p:
#     print(f"{get_current_project()=}")
#     print(f"{p.opened=}")
# print(f"{get_current_project()=}")
# print(f"{p.opened=}")

# projects = get_projects()
# for p in projects:
#     print(p.desc())
#     print(p.opened)
