# Delete Project tasks

If you want to clean up storage used by all tasks within a specific project, you
can loop through them and remove their remote files.This will remove both input
and output files associated with each task.

This is useful when a project is complete or you're starting fresh and want to
free up space.

```python
import inductiva

project = inductiva.projects.Project("my-project")

for task in project.get_tasks():
    task.remove_remote_files()
```
