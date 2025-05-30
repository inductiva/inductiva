# do not import all endpoints into this module because that uses a lot of memory and stack frames
# if you need the ability to import all endpoints from this module, import them with
# from inductiva.client.apis.tag_to_api import tag_to_api

import enum


class TagValues(str, enum.Enum):
    COMPUTE = "compute"
    EVENTS = "events"
    PROJECTS = "projects"
    SIMULATORS = "simulators"
    STORAGE = "storage"
    TASKS = "tasks"
    USERS = "users"
    VERSION = "version"
