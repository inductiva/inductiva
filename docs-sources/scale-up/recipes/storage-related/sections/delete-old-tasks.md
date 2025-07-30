# Delete old tasks

As your project grows, older tasks may accumulate and continue occupying valuable
storage, often without serving any ongoing purpose. If you no longer need the
outputs of these older tasks, you can safely delete their storage to free up
space and reduce costs.

```python
import inductiva
import datetime

# Timezone
timezone = datetime.timezone.utc

# year, month, day, hours, minutes, second
cutoff_datetime = datetime.datetime(2024,7,1,12,34,0,tzinfo=timezone)

# Iterate over all tasks
all_tasks = inductiva.tasks.get_all()

# OR
# Iterate over all tasks from a project
# project = inductiva.projects.Project("my-project")
# all_tasks = project.get_tasks()

for task in all_tasks:
    task_date = task.info.input_submit_time
    if task_date is not None and task_date < cutoff_datetime:
        task.remove_remote_files()
```

```{banner_small}
:origin: recipes_delete_old_tasks
```