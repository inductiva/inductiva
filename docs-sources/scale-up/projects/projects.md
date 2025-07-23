# Projects

Projects are the primary organizational unit in the Inductiva API. Every simulation task belongs to a project, providing a logical way to group related computational work.

## Quick Start

The `Project` class is the core organizational unit in the Inductiva API, providing a simple way to group and manage related simulation tasks. A typical workflow involves three simple steps:
1. **Create** or reference a project for organizing your work
2. **Execute** simulation tasks within the project context
3. **Monitor** and manage all tasks belonging to the project

```python
import inductiva

# Step 1: Create - Create or get reference to a project
project = inductiva.projects.Project("wave_modeling")

# Step 2: Execute - Run simulations and assign to project
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

simulator = inductiva.simulators.XBeach()

# Add the task to the project during its creation
task = simulator.run(input_dir="data/", on=machine_group, project="wave_modeling")

# Or alternatively assign the task after its been created to the project
task = simulator.run(input_dir="data/", on=machine_group)
project.add_task(task)

# Step 3: Monitor - List and manage project tasks
tasks = project.get_tasks()
print(f"Project has {len(tasks)} tasks")

machine_group.terminate()

```

## `Project` Class

A `Project` object serves as a _folder_ for organizing related simulation tasks. Projects are uniquely identified by their name and provide a logical grouping mechanism for computational work.

### Key Features

- **Unique Indentity**: Projects are uniquely identified by their **name** - creating a project with an existing name returns a reference to that project:

```python
import inductiva

# Create or get reference to existing project
project1 = inductiva.projects.Project("wave_modeling")
project2 = inductiva.projects.Project("wave_modeling")

# Different objects, same project
print(project1 is project2)  # False - different objects
print(project1 == project2)  # True - same project
```

- **Default Project**: Every user automatically gets a default project when their account is created. Tasks submitted without specifying a project are automatically assigned to this project.

```python
# These tasks go to the default project
task = simulator.run(input_dir="data/", on=machine_group)
```

### Key Methods

| Method | Purpose | When to Use |
|--------|---------|-------------|
| `Project(name)` | Create or reference a project | Initialize project object |
| `add_task(task)` | Add a task to the project | Explicitly assign existing tasks |
| `get_tasks()` | List all tasks in the project | Monitor project progress |
| `download_outputs()` | Download outputs of all tasks in the project | Analyze project results |

````{eval-rst}
.. seealso::
   For complete API documentation including all parameters, methods, and configuration options, see the `Projects <https://inductiva.ai/guides/api-functions/api/inductiva.projects>`_ class documentation
````


## Create a Project

Create a project by instantiating the `Project` class with a descriptive name:

```python
import inductiva

climate_project = inductiva.projects.Project("wave_modeling")
print(climate_project)
# Output: Project 'wave_modeling' created at 2025-07-22 11:48.
#         Total number of tasks: 0

#         Tasks by status:

#         Estimated total computation cost: 0 US$
```

Or during the task submission process using gthe `run` method:

```python
import inductiva

task = simulator.run(
    input_dir="data/",
    on=machine_group,
    project="wave_modeling"  # if a project with this name does not exist, it will be created
)
```
