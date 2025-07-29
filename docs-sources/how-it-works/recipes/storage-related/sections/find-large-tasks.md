# 🔍 Find Large projects or tasks

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

Now that you have a list of large tasks you can [download and delete](download-and-delete) them.

```{banner_small}
:origin: recipes_find_large_tasks
```