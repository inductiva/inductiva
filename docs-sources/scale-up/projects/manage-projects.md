# Manage Projects

Once you've created projects, you need to manage tasks within them. This guide covers all the essential project management operations and organizational strategies.

| Method | Purpose | Example Usage |
|--------|---------|---------------|
| `simulator.run(project="name")` | Add task during creation | `task = simulator.run(input_dir="data/", on=machine_group, project="coastal_study")` |
| `project.add_task(task)` | Add existing task to project | `coastal_project.add_task(existing_task)` |
| `project.get_tasks()` | List all tasks in project | `tasks = project.get_tasks()` |
| `project.describe()` | Get project summary | `info = project.describe()` |
| `inductiva.projects.get_projects()` | List all projects | `all_projects = inductiva.projects.get_projects()` |

## Adding Tasks to Projects

### Method 1: During Task Creation

The most straightforward approach is to specify the project when submitting a simulation:

```python
import inductiva

# Set up simulator and resources
simulator = inductiva.simulators.XBeach()
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Submit task directly to a specific project
task = simulator.run(
    input_dir="wave_data/",
    on=machine_group,
    project="wave_modeling"  # Task goes directly to project
)

print(f"Task {task.id} submitted to project: {task.get_info().project}")
```

> ⚠️ **Important**: If the project name you specify int the `project` argument does not exist yet, a **_new project with that name will be created_** 

````{eval-rst}
.. ⚠️Important::
   After completing the simulations, remember to release your computational resources 
   to avoid unnecessary charges!
````

### Method 2: After Task Creation

You can also add existing tasks to projects using the `add_task` method:

```python
import inductiva

# Create task without specifying project (goes to default)
task = simulator.run(input_dir="data/", on=machine_group)
print(f"Task initially in project: {task.get_info().project}")  # "default"

# Move task to specific project
wave_modeling = inductiva.projects.Project("wave_modeling")
erosion_project.add_task(task)  # moves to the existing project

# Verify the move
updated_info = task.get_info()
print(f"Task now in project: {updated_info.project}")  # "wave_modeling"
```

## Listing Tasks in Projects

Get all tasks belonging to a project:

```python
import inductiva

project = inductiva.projects.Project("wave_modeling")
tasks = project.get_tasks()
```

Alternatively, you can use the `inductiva projects list` CLI command:

```bash
$ inductiva projects list

 NAME           NR_TASKS
 wave_modeling     2
 ...
 default   241     # default project
```

## Full Example

```python
import inductiva

# Set up computational resources
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Download example input data
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "xbeach-input-example.zip", unzip=True)

# Create project for organized workflow
xbeach_project = inductiva.projects.Project("wave_modeling")
print(f"Working with project: {xbeach_project.name}")

# Initialize simulator
simulator = inductiva.simulators.XBeach()

# Method 1: Submit tasks directly to project
task_1 = simulator.run(
    input_dir=input_dir,
    on=machine_group,
    project="xbeach_project"
)

# Method 2: Create task first, then assign to project
task_2 = simulator.run(input_dir=input_dir, on=machine_group)
xbeach_project.add_task(task_2)

# List all tasks in the project
project_tasks = xbeach_project.get_tasks()
print(f"\nProject has {len(project_tasks)} tasks:")

for task in project_tasks:
    print(f"  • {task.id}: {task.get_status()}")
    print(f"    {task.print_summary()}")

# Clean up resources
machine_group.terminate()
```