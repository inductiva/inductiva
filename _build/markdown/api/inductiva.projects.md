# projects

This module provides functionality for managing projects and their associated
tasks within the Inductiva platform. A project serves as a container for
grouping related tasks, enabling better organization and management of
computational workflows.

Classes:
: - Project: Represents a project that groups related tasks together.

Functions:
: - get_projects(): Retrieves all projects associated with the current user.

Key Features:
: - Create and manage projects on the backend.
  - Add tasks to projects and retrieve the most recent tasks or filter by
    status.
  - Monitor the status of all tasks in a project and wait for their
  <br/>
  completion.
  - Download outputs for all tasks in a project.
  - Estimate the total computation cost of a project based on its tasks.

Example Usage:

> ```python
> import inductiva

> # Create a new project or load an existing one
> project = inductiva.projects.Project("my_project")

> # Select a simulator
> fvcom = inductiva.simulators.FVCOM()
> # Run and add tasks to the project
> for i in range(10):
>     task = fvcom.run(...)
>     project.add_task(task)

> # Monitor task completion
> project.wait()
> # Download outputs for all tasks
> project.download_outputs()
> # Print project details
> print(project)
> ```

### *class* Project(name: str)

Bases: `object`

Projects management class.

Groups related tasks together under a single project.

### Example

```python
project = inductiva.projects.Project("my project")

task_1 = simulator.run(...)
project.add_task(task_1)

task_2 = simulator.run(...)
project.add_task(task_2)
```

#### \_\_init_\_(name: str)

Initialize the Project instance.

* **Parameters:**
  **name** (*str*) – The name of the project.

#### add_task(task: [Task](inductiva.tasks.md#inductiva.tasks.task.Task))

Adds a task to the project.

* **Parameters:**
  **task** – The task to add to the project.

#### *property* created_at *: str*

Returns the creation date and time of the project.

#### delete()

Delete a project on the backend.

This method does not delete the project tasks, only the project itself.
The tasks will be moved to the “default” project.

#### download_outputs(output_dir: str | None = None)

Downloads all the outputs for all the tasks in the project.

All task outputs will be organized within the specified output_dir.
If output_dir is not provided, outputs will be saved to a default
location under inductiva_output/<project_name>/<task_id>/.
Otherwise, they will be stored in <output_dir>/<task_id>/.

* **Parameters:**
  **output_dir** (*str* *,* *optional*) – The base directory where project outputs
  will be downloaded.

#### *property* estimated_computation_cost *: float*

Returns the estimated project cost.

Computed as the sum of the estimated computation cost of each task.

#### get_tasks(last_n: int = -1, status: str | None = None) → List[[Task](inductiva.tasks.md#inductiva.tasks.task.Task)]

Get the the tasks of this project.

Optionally, those can be filtered by task status.

* **Parameters:**
  * **last_n** (*int*) – The number of tasks with repect to the submission
    time to fetch. If last_n<=0 we fetch all tasks submitted
    to the project.
  * **status** – Status of the tasks to get. If None, tasks with any
    status will be returned.

#### *property* id *: str*

Returns the unique ID of the project.

#### *property* name *: str*

Returns the name of the project.

#### *property* num_tasks *: int*

Returns the number of tasks in the project.

#### *property* task_by_status *: dict*

Returns a dictionary with the number of tasks by status.
The keys are the status codes and the values are the number of tasks
with that status.

#### *property* total_estimated_cost *: float*

#### *property* total_task_orchestration_fee *: float*

#### wait()

Wait for all the tasks in a project to complete.

### get_projects() → List[[Project](#inductiva.projects.project.Project)]

Gets all the user’s projects.

::docsbannersmall
::
