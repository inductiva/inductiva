# Tasks

Central to the API's functionality is the concept of a `Task`,
which is an  abstraction that encapsulates all the information
about the computational workload you defined as the user.
Once you submit a simulation to the  API, a `Task` is generated.
This allows for real-time updates on simulation status, including
monitoring its progress and retrieving its outputs.

In this reference, you will learn about the entire lifecycle of a `Task`,
including its creation, various operational states, and termination. This information
is crucial to understand how the Inductiva API manages simulations.

## Task Creation

A `Task` is created when you submit a simulation via the API by invoking the `run`
method on a simulator object. Each call to this method, even with identical arguments,
generates a unique `Task`, leading to separate executions of the same simulation.

Here's an example of how to create a `Task` by submitting a SpliSHSPlasH simulation
to the API:

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

Note that a subsequent call to `splishsplash_simulator.run()` with the same
input arguments would create a new, distinct task:

```python
task2 = splishsplash_simulator.run(input_dir="splishsplash-input-example",
                                  sim_config_filename="config.json",
                                  on=machine_group)
print (task2.id)
# Example output: k9muu1vq1fc6m2oyxm0n3n8y0

machine_group.terminate()
```

Each `Task` is identified by a unique alphanumeric identifier. While you can
instantiate multiple `Task` objects pointing to the same task using this identifier,
this does not duplicate the task on the API. These objects refer to the same underlying
task. This mechanism is useful when you want to recreate a `Task` object to
retrieve information about tasks you created in previous sessions:

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

## Task Execution

Tasks run **asynchronously** by default, enabling you to continue interacting with
the API without interrupting your code in the background, even as tasks run or are
queued for execution. This asynchronous model also allows for batch submissions within
a single session without waiting for each to complete, with the API efficiently
distributing them across available computational resources for parallel processing.

You can configure how to wait for the task to end by either using the `wait` method or
the `sync_context` manager.

**`wait`**

You can turn a task into a blocking call by using the `wait` method.
This method will block the call, allowing your script to wait until the task comes
to a terminal status without affecting the remote simulation if it were locally
interrupted (ex. the script is abruptly terminated).

This is useful when you want to wait for the completion of a task before proceeding
with some other operation (ex. downloading the output files):

```python
# Block the call until the simulation completes
task.wait()  # <- The remote simulation WILL NOT DIE if the local session is
             #    interrupted while waiting for the wait() call to return
```

**`sync_context`**

Alternatively, you can use the `sync_context` manager to ensure the remote simulation
terminates if your local session is interrupted while waiting for the `wait` call
to return (e.g., through a Ctrl+C command):

```python
# Block the call until the simulation completes
with task.sync_context():
    task.wait()  # <- The remote simulation WILL DIE if the local session is
                 #     interrupted while waiting for the wait() call to return
```

(task-lifecyle)=
## Task Lifecycle

The status of a task changes as it progresses through its lifecycle. Understanding
these states is crucial if you want to track a task's progress through the API and
manage your simulations effectively.

The following diagram shows the path a task may take through the API, and identifies
the relevant states and possible state transitions:


<div align="center">
   <img src="../_static/task_state.svg" alt="Task state diagram">
   <figcaption align = "center"><b>State diagram of a Task</b></figcaption>
</div>

---

Below a succinct description of each state, including the actions that
lead to a state transition:

|  	|  	|
|---	|---	|
| `PENDING INPUT` 	| After you make your initial simulation request, the newly created task remains in this pending state, awaiting all necessary input files. You can terminate or kill the task, moving it to `KILLED`, or it may become `ZOMBIE` if its assigned machine group is terminated. Once you've uploaded all the input files, the task progresses to `SUBMITTED`. 	|
| `SUBMITTED` 	| The task is queued and ready to be picked up by an executor. You can cancel the task, moving it to `KILLED`, or it can become a `ZOMBIE` under the same conditions as in `PENDING INPUT`. When an executor picks it up, the task moves to `STARTED`. 	|
| `STARTED` 	| Simulation has started. Upon successful completion, the task transitions to `SUCCESS`. If it fails due to simulation issues, the task transitions to `FAILED`,whereas any failure due to executor problems moves the task to `EXECUTOR FAILED`. You can send a request to kill the task, moving it to `PENDING KILLED`. 	|
| `PENDING KILLED` 	| Your request to terminate a running task has been received by the API and is awaiting execution. 	|
| `KILLED` 	| The task has been successfully terminated upon your request. 	|
| `ZOMBIE` 	| The progression of the non-started task abruptly stops due to the shutdown of the computational resources where the task was running, typically when a machine group is user-terminated. 	|
| `SPOT PREEMPTED` 	| Spot instances running the task were terminated. In this case, the task is re-queued or `SUBMITTED` until new resources with the same original machine group become available. 	|
| `FAILED` 	| Simulator errors prevent completion, usually due to incorrect input configurations or an internal error within the simulator itself. 	|
| `EXECUTOR TERMINATED` 	| The executor was terminated due to internal reasons. Similar to `SPOT PREEMPTED`, the task is requeued or `SUBMITTED` for execution. 	|
| `EXECUTOR TERMINATED BY USER` 	| Similar to `KILLED`, the task's executor was terminated by you. 	|
| `EXECUTOR TERMINATED TTL EXCEEDED` 	| The executor running the task was terminated due to exceeding its maximum time to live.   |
| `EXECUTOR FAILED` 	| Executor errors prevent completion, such as low disk space. 	|
| `TTL EXCEEDED` 	| The task exceeded its configured time to live, *i.e.*, a limit to the task's computation time when running in a shared resource.  |
