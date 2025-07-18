# Projects

Projects are the primary way to organize `Tasks` in the Inductiva API.
A project is a collection of tasks that run under a common umbrella, identified
by its name. Any submitted task will always belong to a project: either a
user-specified one or the default project. This allows users to organize their
tasks in a way that makes sense to them.
The following sections will explain how to create and interact with projects.

## Creating a Project

A project is created whenever an instance of the `Project` class, from the
`projects` module, is created. **Projects are uniquely identified by their name**.

When instantiating a new project, the user needs to specify the project's name -
if one with the same name already exists, a reference to that project
will be returned, otherwise, a new one is created. The following snippet
creates an empty project called `my_project` and prints a description of it:

```python
>>> import inductiva
>>> project = inductiva.projects.Project("my_project")
>>> print(project.describe())
Project 'my_project' with 0 tasks (id=35ea9cf0-e8ce-47cf-b6cf-6a3eb4994d98)

# At this point, any new Project instantiated with the same name will
# refer to the same project, even though the object is different:
>>> project2 = inductiva.projects.Project("my_project")
>>> print(project2 is project) # test if the objects are the same
False
>>> print(project2 == project) # test if the projects are the same
True
```

**NOTE:** A 'default' project is automatically created when a new user account is added.
This project is used when no project is explicitly used in the task submission.

## Listing existing projects

To retrieve a list of all existing projects, use the `get_projects` function
from the `projects` module:

```python
>>> import inductiva
>>> projects = inductiva.projects.get_projects()
[<inductiva.projects.project.Project at 0x123f22110>,
 ...
 <inductiva.projects.project.Project at 0x123f00190>]
```

Alternatively, you can use the `inductiva projects list` CLI command to list
all existing projects:

```bash
$ inductiva projects list

 NAME           NR_TASKS
 my_project     0
 ...
 userab1cdef2   241     # default project
```

## Adding a Task to an existing Project

To add a task to an existing project, e.g. "demo", the user can simply pass
the name of the project as an argument to the ```simulator.run``` method.

```python
task = simulator.run(input_dir=input_dir,
                     on=machine_group,
                     project="demo")

```

Alternatively, if you have a `Task` and a `Project` objects, you can call the `add_task`
method explicitely:

```python

my_demo_project.add_task(task_1)

```


### Full Example

The following snippet demonstrates how to add a task to a new project
and another one to the default project.

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# get example input data
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "xbeach-input-example.zip", unzip=True)

project = inductiva.projects.Project("my_xbeach_project")

simulator = inductiva.simulators.XBeach()

task1 = simulator.run(input_dir=input_dir,
                      on=machine_group)

# add a task to the "my_xbeach_project" project
project.add_task(task_1)

# task2 will be added to the default project
task2 = simulator.run(input_dir=input_dir,
                      on=machine_group)

print(task1.get_info().project) # "my_xbeach_project"
print(task2.get_info().project) # "default"

machine_group.terminate()
```

### [Deprecated] Using a Context Manager

Until inductiva version 0.15.9 Projects needed to be "opened" and "closed"
and there was a syntax using Python context managers to call these operators
under the hood.
These mechanisms were deprectated with version 0.16.0, in favour of the simpler
alternatives documented above. 

**DEPRECATED**
```python

with inductiva.projects.Project("my_xbeach_project", append=True) as project:
    simulator = inductiva.simulators.XBeach()

    # add a task to the "my_xbeach_project" project
    task1 = simulator.run(input_dir=input_dir,
                          on=machine_group)
```

## Listing Tasks in a Project

Each project has a list of tasks that belong to it. To list all tasks in a given
object, use the `get_tasks` method of the `Project` class:

```python
>>> import inductiva
>>> project = inductiva.projects.Project("my_xbeach_project")
>>>> tasks = project.get_tasks()
[<inductiva.tasks.task.Task at 0x123f22110>,
 ...
 <inductiva.tasks.task.Task at 0x123f00190>]
```

Alternatively, one can use the `inductiva tasks list` CLI command to list all
tasks of a project using the `-p/--project` filter argument:

```bash
$ inductiva tasks list -p my_xbeach_project
Showing tasks for project: my_xbeach_project.

 ID                          SIMULATOR   STATUS     SUBMITTED          STARTED            COMPUTATION TIME     RESOURCE TYPE
 2qnbmu6jxnf4zv8ma0c19ujhe   xbeach      success    24 May, 15:31:37   24 May, 15:31:45   0:00:01              GCP c2-standard-4
 [...]
 q3k3ad1etqdwwfw31die00orw   xbeach      success    24 May, 13:11:21   24 May, 13:11:21   0:00:06              GCP c2-standard-4

```
