# Delete failed tasks

Over time, failed tasks can accumulate and take up unnecessary storage space. To
free up resources, you might want to remove the storage associated with these
failed runs.

```python
import inductiva

# Iterate over all tasks
all_tasks = inductiva.tasks.get_all(status="failed")

# Iterate over all tasks from a project
# project = inductiva.projects.Project("my-project")
# all_tasks = project.get_tasks(status="failed")

for task in all_tasks:
    task.remove_remote_files()
```
> **Note**: You can query tasks by many different status. Learn more about a task lifecycle [here](https://inductiva.ai/guides/how-it-works/intro/tasks#task-lifecycle).
