from inductiva.client.paths.tasks_task_id_metadata.get import ApiForget
from inductiva.client.paths.tasks_task_id_metadata.put import ApiForput
from inductiva.client.paths.tasks_task_id_metadata.patch import ApiForpatch


class TasksTaskIdMetadata(
        ApiForget,
        ApiForput,
        ApiForpatch,
):
    pass
