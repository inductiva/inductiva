# Delete all tasks from a project.

If you want to clean up storage used by all tasks within a specific project, you
can loop through them and remove their remote files. This is useful when a
project is complete or you're starting fresh and want to free up space.

```python
import inductiva

project = inductiva.projects.Project("my-project")

for task in project.get_tasks():
    task.remove_remote_files()
```

# Identify Projects or Tasks That Require Large Storage Volumes.

It’s not always obvious which projects or tasks are consuming the most storage
space. To help you pinpoint what’s taking up significant resources, here’s a
simple method to identify tasks or projects that exceed a given storage
threshold (e.g., X GB).

```python
import inductiva

# Iterate over all tasks, regardless of project
all_tasks = inductiva.tasks.get_all()

# OR
# Iterate over all tasks from a project
# project = inductiva.projects.Project("my-project")
# all_tasks = project.get_tasks()

# Threshold in GB
threshold = 2.1

big_tasks = []

for task in all_tasks:
    # If the task has no storage we set it to 0
    task_storage = task.info.data_metrics.output_zipped_size_bytes.value or 0
    
    gb_storage = task_storage / (1024 ** 3)
    
    if gb_storage>= threshold:
        print(f"The task {task.id} has an output of {gb_storage:.2f} GB.")
        big_tasks.append(task)
```

This script populates the `big_tasks` list with all tasks whose outputs exceed
2.1 GB. You can adjust the `threshold` value to whatever limit best fits your needs.

You can now do whatever you want to the list of flagged tasks.

## Download and delete the storage from a list of tasks.

You now have a list of tasks that you consider that are too big to sit in your storage
eating space and melting your credits. Let's do a simple loop where you download
the output files and after that you delete them from your remote storage.

```python
import inductiva

for task in big_tasks:
    task.download_outputs()
    task.remove_remote_files()
```

It's this simple.

> **Note**: Keep in mind that downloading the files has a cost (learn more about costs [here](https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost)). If you can,
download only the needed files instead of the whole storage. (learn more how you
can do that [here](download-file-from-project).)

# How can I delete the storage from all failed tasks?

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

# How can I delete the storage from all old tasks?

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