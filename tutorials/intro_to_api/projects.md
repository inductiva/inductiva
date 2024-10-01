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

**NOTE:** A default project is automatically created when a new user account is added.
Typically, the default project is named after the user's username with some random
characters appended to it. This project is used when no project is explicitly used
in the task submission and ensures that no orphaned tasks are created.

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

## Adding a Task to a Project

Any task that is submitted to the Inductiva API must belong to a project.
When the user is not explicit about the project to which the task is to be added,
the default one is used. This ensures that no orphaned tasks are created.
To add a task to a project other than the default one, the user needs to create
a new project (or reference an existing one) and open it for task submission. However,
the project can only be opened for task submission if it gets instantiated with
an explicit `append=True` argument in the constructor or, otherwise, a RunTime
error will be raised. This mechanism is meant to work as a security measure that
ensures that the user is explicitly aware that new tasks are being appended to the
project, which is particularly useful to ensure the consistency of the project's
content.

```python
>>> import inductiva
>>> project = inductiva.projects.Project("demo", append=False)
>>> project.open() # <-- An exception will be raised because the project is appendable
...
RuntimeError: Trying to open a project with `append=False`.
A Project can only be opened when instantiated with the `append=True` option.
```

Whenever the `run` method of a `Simulator` object is called, the resulting task
will be created under the project that is currently open, even though the project
is not explicitly specified as an argument to the `run` method. It will be inferred
from the context where the `run` method is being called.

The client library offers two mechanisms to manage how tasks are added to a
project: by manually opening and closing a `project` or via a context manager.

### Explicit management

Explicit management requires that the user opens and closes the project manually
using the `open` and `close` methods of the `Project` class. The project needs to
be instantiated, opened for task submission and then closed to prevent further tasks
from being appended to it. Any task submitted between the opening and closing of
the project will be added to it. The following snippet demonstrates how to add a
task to a project using explicit management:

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup('c2-standard-4')
machine_group.start()

# get example input data
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "xbeach-input-example.zip", unzip=True)

project = inductiva.projects.Project("my_xbeach_project", append=True)

project.open() # <-- open the project for task submission

simulator = inductiva.simulators.XBeach()

# add a task to the "my_xbeach_project" project
task1 = simulator.run(input_dir=input_dir,
                      sim_config_filename="params.txt",
                      on=machine_group)

project.close() # <-- close the project

# task2 will be added to the default project
task2 = simulator.run(input_dir=input_dir,
                      sim_config_filename="params.txt",
                      on=machine_group)

print(task1.get_info().project) # "my_xbeach_project"
print(task2.get_info().project) # "userab1cdef2" (default project)

machine_group.terminate()
```

### Using a Context Manager

Alternatively, the project can be managed using a context manager that automatically
opens and closes the project for task submission, and helps clarify the scope of
the project. The manager ensures that the project is opened for task submission
by calling the `open` and `close` methods when entering and exiting the context
of the `with` block, respectively. Any task that is created inside the `with` block
will be appended to the project managed in that context; tasks submitted outside
will be added to the default one:

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup('c2-standard-4')
machine_group.start()

# get example input data
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "xbeach-input-example.zip", unzip=True)

with inductiva.projects.Project("my_xbeach_project", append=True) as project:
    simulator = inductiva.simulators.XBeach()

    # add a task to the "my_xbeach_project" project
    task1 = simulator.run(input_dir=input_dir,
                          sim_config_filename="params.txt",
                          on=machine_group)

# task2 will be added to the default project
task2 = simulator.run(input_dir=input_dir,
                      sim_config_filename="params.txt",
                      on=machine_group)

print(task1.get_info().project) # "my_xbeach_project"
print(task2.get_info().project) # "userab1cdef2" (default project)

machine_group.terminate()
```

At any moment, the user can query what project is currently **open** for task submission
by calling the `get_current_project` function from the `projects` module:

```python
>>> inductiva.projects.get_current_project()
None
>>> with inductiva.projects.Project("my_xbeach_project", append=True):
...     inductiva.projects.get_current_project().name
my_xbeach_project
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

## Thread awareness of Projects

Projects are thread-aware. Multiple threads can define and use different
projects simultaneously. A project opened in one thread will not affect the
project opened in another thread. If a thread does not open a project, the default
one will be used, even though the main thread might have opened a project other
than the default one.

At any point, the user can get a reference to the currently active project in
the calling thread using the `inductiva.projects.get_current_project()` function.
This function will return the project that is currently open for task submission,
or `None` if the default project is being used. The following snippet demonstrates
how to use this function:

```python
>>> import inductiva
>>> # so far, no project is open, hence "None" is returned
>>> inductiva.projects.get_current_project() is None
True
>>> # open the "demo" project and get a reference to it
>>> with inductiva.projects.Project('demo', append=True) as p:
...     print(inductiva.projects.get_current_project().name)
...     print(inductiva.projects.get_current_project().num_tasks)
...     print(inductiva.projects.get_current_project() is p)
demo
4
True
>>> # the "demo" project is closed, hence "None" is returned again
>>> inductiva.projects.get_current_project() is None
True
```

To illustrate the thread-awareness of projects, consider the following example.
The main thread opens a project, and a new thread is created in the context of the
opened project. The new thread retrieves the current project and prints its object
just to find out that the default project is being used therein, even though
the main thread has opened a project other than the default one:

```python
import inductiva
import threading

def run():
    print(f"Thread 1: {inductiva.projects.get_current_project()=}")

with inductiva.projects.Project("demo", append=True):
    print(f"Main thread: {inductiva.projects.get_current_project()=}")
    thread = threading.Thread(target=run)
    thread.start()
    thread.join()

# would print:
# Main thread: inductiva.projects.get_current_project()=<inductiva.projects.project.Project object at 0x1242be910>
# Thread 1: inductiva.projects.get_current_project()=None
```
