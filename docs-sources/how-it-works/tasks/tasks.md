# Tasks

A `Task` is automatically created whenever you submit a simulation to the API by calling the `run` method on a simulator object. Each call to this method generates a unique `Task`, even if the arguments remain identical. This ensures separate executions for each submission of the same simulation.

## Submitting Your First Task

Here's an example of how to create a `Task` by submitting a [SpliSHSPlasH](https://inductiva.ai/guides/splishsplash) simulation to the API:

```python
# Running the SplishSplash simulator example
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Download example input files
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "splishsplash-input-example.zip", unzip=True)

# Instantiate the simulator
splishsplash_simulator = inductiva.simulators.SplishSplash()

# Run the SPlisHSPlasH simulation with the required .json config file
task = splishsplash_simulator.run(input_dir="splishsplash-input-example",
                                  sim_config_filename="config.json",
                                  on=machine_group)
print(task.id)
# Example output: i4ir3kvv62odsfrhko4y8w2an
```

Note that a subsequent call to `splishsplash_simulator.run()` with the same input arguments would create a new, distinct task:

```python
task2 = splishsplash_simulator.run(input_dir="splishsplash-input-example",
                                  sim_config_filename="config.json",
                                  on=machine_group)
print (task2.id)
# Example output: k9muu1vq1fc6m2oyxm0n3n8y0

machine_group.terminate()
```

When a task is submitted, the API provides immediate feedback in the terminal about its position in the queue and helpful CLI commands for monitoring:

```sh
# Number of tasks ahead in the queue: 3.
# · Consider tracking the status of the task via CLI:
#     inductiva tasks list --id <task_id>
# · Or, tracking the logs of the task via CLI:
#     inductiva logs <task_id>
# · Or, track the task files in real time with:
#     inductiva tasks list-files <task_id>
# · You can also get more information about the task via the CLI command:
#     inductiva tasks info <task_id>
```

You can programmatically check the queue position and task status:

```python
# Check queue position
position = task.get_position_in_queue()
if position is not None:
    print(f"Tasks ahead in queue: {position}")

# Check current status
print(f"Current status: {task.status}")
```

## Unique Task Identification

Every `Task` is assigned a unique alphanumeric identifier upon creation. This identifier ensures that each task is distinct and can be easily referenced.

While it is possible to instantiate multiple `Task` objects using the same identifier, doing so does not create duplicate tasks on the API. Instead, all such objects point to the same underlying task, allowing you to access and manage it across different sessions.

This mechanism is useful when you want to recreate a `Task` object to retrieve information about tasks you created in previous sessions:

```python
>>> import inductiva
>>>
>>> task1 = inductiva.tasks.Task("i4ir3kvv62odsfrhko4y8w2an")
>>> task2 = inductiva.tasks.Task("i4ir3kvv62odsfrhko4y8w2an")
>>> print(id(task1))
4410160112
>>> print(id(task2))
4389863104
>>> print(task1.id)
i4ir3kvv62odsfrhko4y8w2an
>>> print(task2.id)
i4ir3kvv62odsfrhko4y8w2an # Outputs the same task ID as task1.id
```

## Task Metadata

Many Inductiva API users generate extensive datasets by running hundreds or thousands of
simulation tasks, each with different input parameters. Previously, tracking these parameters
required maintaining external mechanisms — either local files or databases — to associate task IDs
with their corresponding parameter sets.

Now, it is possible to store these parameters as `task metadata` directly on the Inductiva platform:

```python
# Run a SpliSHSPlasH simulation
task = splishsplash_simulator.run(input_dir="splishsplash-input-example",
                                  sim_config_filename="config.json",
                                  on=machine_group)

# Store simulation parameters as metadata
task.set_metadata({
    "fluid_viscosity": "0.01",
    "particle_radius": "0.025",
    "simulation_time": "5.0",
    "run_number": "2",
})
```

> Note: both metadata keys and values must be **strings**. Any numeric values, booleans, or other data types should be converted to strings before being stored as metadata.

And the metadata can later be retrieved by doing:

```python
metadata = task.get_metadata()
```

You can also manipulate the metadata in the Task's page on the Web Console:

<div align="center">
   <img src="../_static/task-metadata.png" alt="Task Metadata">
   <figcaption align="center"><b>Task Metadata</b></figcaption>
</div>