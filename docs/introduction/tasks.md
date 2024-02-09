# Tasks

At the heart of the API lies the concept of a `Task`.
A task is an abstraction that encapsulates all the information about the
computational workload defined by the user. A task is generated when a
simulation is submitted to the API and exists until the user deletes it.
With a task at hand, the user can manage the simulation, retrieve its status,
download its output files, and more.

In this section, we will cover the following topics:
- creating a task, including retrieving tasks from previous sessions;
- managing a task;
- monitoring a task;
- using the task to retrieve the results of a simulation.


## Task Creation

A task is created when a simulation is submitted to the API through a call to the
`run` method of a simulator object. _N_ calls to that method using the same
arguments will produce _N_ tasks, all referring to distinct executions of the
same simulator to perform the same simulation. In the following snippet, we show
how to create a task by submitting a SpliSHSPlasH simulation to the API.

```python
# Example of how to run the SplishSplash simulator
import inductiva

# get some example input files
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "splishsplash-input-example.zip", unzip=True)

# Instantiate the simulator
splishsplash_simulator = inductiva.simulators.SplishSplash()

# SPlisHSPlasH requires a .json file as the main config file.
task = splishsplash_simulator.run(input_dir="splishsplash-input-example",
                                  sim_config_filename="config.json")
print(task.id)  # outputs i4ir3kvv62odsfrhko4y8w2an (example id)

# Note that the following call to splishsplash_simulator.run() using the same
# input arguments would create a new and distinct task
# task2 = splishsplash_simulator.run(input_dir="splishsplash-example",
#                                  sim_config_filename="config.json")
# print (task2.id)  # outputs k9muu1vq1fc6m2oyxm0n3n8y0 (example id)

```

A task is uniquely identified by an alphanumeric identifier which is unique across
all tasks submitted to the API by the user. However, a user can create multiple
`Task` objects that point to the same task by instantiating the `Task` class with
the same identifier, without this resulting in the creation of new tasks on the
API; they will be different Python objects, but they will refer to the same task
on the API. This mechanism is useful when the user wants to recreate a `Task`
object to retrieve information about tasks created in previous sessions.

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
i4ir3kvv62odsfrhko4y8w2an  # the same as task1.id
```

### Synchronous tasks

Tasks are asynchronous by default. This means that the user can continue to
interact with the API while the task is running or queued up for execution.
This behavior allows the user to submit batches of tasks in a single session
knowing that the API will distribute them across the available computational
resources so that they can run in parallel.

However, the user can turn a task into a blocking call by using the `wait` method.
This method will block the call until the task comes to a terminal status, but
without interrupting the remote simulation if the local session is interrupted
(ex. the script is abruptly terminated). This is useful when the user wants to
wait for the completion of a task before proceeding with some other operation
(ex. downloading the output files).

```python
# Block the call until the simulation completes
task.wait()  # <- The remote simulation WILL NOT DIE if the local session is
             #    interrupted while waiting for the wait() call to return
```

Alternatively, the user can turn the simulation into a blocking call by using
the `sync_context`context manager. This context manager ensures that the remote
simulation is indeed killed if the local session is interrupted while waiting
for the `wait` call to return.

```python
# Block the call until the simulation completes
with task.sync_context():
    task.wait()  # <-- The remote simulation WILL DIE if the local session is
                 #     interrupted while waiting for the wait() call to return
                 #     Ex. the user presses Ctrl+C while waiting for the call to
                 #     return
```

This distinction is relevant because, at any moment, the user might wonder
whether interrupting the local session will also kill the remote simulation.
A task will only be killed if the user explicitly requests it or if the session
is terminated while waiting for a `wait` call to return within the context of a
`sync_context()` call.


### Task Lifecycle


The status of a task changes as it progresses through its lifecycle.
So far, we have given an overview of task creation and, to some extent,
task killing but yet, a myriad of intermediate states exist, which are
important for the user to track the progress of the request through the API.
Thus, tasks are characterized by their state.

The following diagram shows the lifecycle of a task, and identifies the relevant states and possible state transitions:


<div align="center">
   <img src="../_static/task_state.svg" alt="Task state diagram">
   <figcaption align = "center"><b>State diagram of a Task</b></figcaption>

</div>


We now provide a succinct description of each state, including the actions that
lead to a state transition:

- `PENDING INPUT`: upon issuing a request for a simulation by the user, the newly
    created task is in this state. It remains in this state until all input files
    are uploaded to the API. The user can kill the task (transitions to `KILLED`)
    or the task can become zombie if the machine group is meanwhile terminated
    (transitions to `ZOMBIE`).
    After all input files are uploaded, the task transitions to `SUBMITTED`.

- `SUBMITTED`: In this state, the task is ready to be picked up by an executor.
    The user can kill the task (transitions to `KILLED`) or the task can become
    zombie if the machine group is meanwhile terminated (transitions to `ZOMBIE`).
    After being picked up by an executor, the task transitions to `STARTED`.

- `STARTED`: The task has been picked from the waiting queue and is being executed.
    Upon successful completion, the task transitions to `SUCCESS`. If it fails
    for some reason pertaining to the simulator, the task transitions to `FAILED`
    whereas if the failure is due to the executor, the task transitions to
    `EXECUTOR FAILED`. The user can send a request to kill the task
    (transitions to `PENDING KILLED`).

- `PENDING KILLED`: The user has requested a running task to be killed and the
    request has been received by the API but has not yet been delivered and
    executed by the executor.

- `KILLED`: The task is killed.

- `ZOMBIE`: The progression of the non-started task was suddenly stopped because
    the computational resources where the task was running where brought to a halt.
    This can happen when the machine group is terminated by the user.

- `SPOT PREEMPTED`: the task was running in spot instances that were terminated.
    In this case, the task is resubmitted for execution whenever new resources
    are available with the same original machine group
    (transitions to `SUBMITTED`);

- `EXECUTOR TERMINATED`: The executor was terminated due to internal reasons.
    Similarly to the `SPOT PREEMPTED` state, the task is resubmitted for
    execution;

- `EXECUTOR TERMINATED BY USER`: The task was running in an executor but was
    killed either by the user or by the system (generally when quota limits are
    reached).

- `FAILED`: The simulation failed to reach completion due to an error in the
    simulation itself (generally due to the bad/miss-configured input files or an
    internal error of the simulator).
    
- `EXECUTOR FAILED`: The simulation failed to reach completion due to an error in
    the executor itself (e.g. no space in disk, etc).

    