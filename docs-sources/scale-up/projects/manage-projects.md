# Manage Projects

Once you've created projects, you need to manage tasks within them. This guide covers all the essential project management operations and organizational strategies.

| Method | Purpose | Example Usage |
|--------|---------|---------------|
| `simulator.run(project="name")` | Add task during creation | `task = simulator.run(input_dir="data/", on=machine_group, project="coastal_study")` |
| `project.add_task(task)` | Add existing task to project | `coastal_project.add_task(existing_task)` |
| `project.get_tasks()` | List all tasks in project | `tasks = project.get_tasks()` |
| `project.wait()` | Wait for all the tasks in a project to complete | `coastal_project.wait()` |
| `project.download_outputs()` | Downloads the outputs for all the tasks in the project | `coastal_project.download_outputs(output_dir="coastal_project_results")` |
| `project.delete()` | Delete the project and moves the tasks to the default project | `coastal_project.delte()` |
| `inductiva.projects.get_projects()` | List all projects | `all_projects = inductiva.projects.get_projects()` |

## Add Tasks to Projects

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

> ⚠️ If the project name you specify int the `project` argument does not exist yet, a **_new project with that name will be created_** 

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

## List Tasks in Projects

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
 coastal_study     5
 default   241     # default project
```

## Wait for a Project

Use the `wait()` method to block execution until all tasks in a project complete:

```python
import inductiva

project = inductiva.projects.Project("wave_modeling")

print("Starting project execution...")
project.wait()  # Blocks until all tasks are done
print("All project tasks completed!")

# Check final status of all tasks
tasks = project.get_tasks()
for task in tasks:
    print(f"Task {task.id}: {task.get_status()}")
```

## Download a Project Outputs

Download all outputs from a project's tasks in one operation:

```python
import inductiva

project = inductiva.projects.Project("wave_modeling")

# Wait for completion first
project.wait()

# Download all outputs to a directory
project.download_outputs(output_dir="wave_modeling_results")
print("All project outputs downloaded to 'wave_modeling_results/'")
```

The downloaded directory structure will organize outputs by task:

```
wave_modeling_results/
├── task_abc123/
│   ├── output_file1.dat
│   └── output_file2.log
└── task_def456/
    ├── output_file1.dat
    └── output_file2.log
```

## Delete a Project

When you no longer need a project, you can delete it using the `delete()` method.

> ⚠️ Deleting a project **moves all its tasks back to the default project** - the tasks themselves are **not deleted**.

```python
import inductiva

# Get the project you want to delete
project = inductiva.projects.Project("old_study")

# Delete the project (tasks move to default project)
project.delete()
print(f"Project '{project.name}' deleted. Tasks moved to 'default' project.")

# Verify tasks moved to default project
default_project = inductiva.projects.Project("default")
default_tasks = default_project.get_tasks()
print(f"Default project now has {len(default_tasks)} tasks")

```

## List All Your Projects

View all your projects:

```python
import inductiva

# Get all projects
all_projects = inductiva.projects.get_projects()

print(f"You have {len(all_projects)} projects:")
for project in all_projects:
    task_count = len(project.get_tasks())
    print(f"  • {project.name}: {task_count} tasks")
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
    project="wave_modeling"
)

# Method 2: Create task first, then assign to project
task_2 = simulator.run(input_dir=input_dir, on=machine_group)
xbeach_project.add_task(task_2)

# List all tasks in the project
project_tasks = xbeach_project.get_tasks()
print(f"\nProject has {len(project_tasks)} tasks:")

# Wait for all tasks to complete
print("\nWaiting for project completion...")
xbeach_project.wait()

# Download all results
print("Downloading project outputs...")
xbeach_project.download_outputs(output_dir="wave_modeling_results")

# Final project summary
print(f"\nProject Summary:")
print(xbeach_project)

# Clean up resources
machine_group.terminate()
```