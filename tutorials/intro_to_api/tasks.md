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

Status | What happened? | Why? | Recommendations to overcome
-- | -- | -- | --
PENDING INPUT | Your simulation is waiting for all necessary input files. | The task can’t start because required input files are missing. If multiple tasks are queued, the first task will process sequentially once the inputs are ready. | Upload the missing input files to proceed or terminate the task if no longer needed. For faster task processing, submit tasks to separate machine groups.
SUBMITTED | Your task is queued and waiting for an executor to start it. | The task is ready but hasn’t yet been picked up by an available machine. If you’ve queued multiple tasks, they will execute one at a time. | Wait for the task to start, or cancel it if no longer needed. To process tasks faster, submit each to a separate machine group.
STARTED | Your simulation has started running. | The task has been picked up by an executor and is now processing. | Wait for the simulation to finish and transition to SUCCESS. If issues arise, check for a FAILED or EXECUTOR FAILED status, or kill the task if no longer needed.
COMPUTATION STARTED | The task's inputs and container image have been downloaded. | The machine is being set up with the appropriate assets to run your task. | The task's computation will now start.
PENDING KILLED | Your request to terminate the task has been received and is awaiting execution. | The system is processing your termination request for the running task. | Wait for the task to be fully terminated.
KILLED | The task has been successfully terminated as per your request. | The task was killed either because it exceeded the set time-to-live (TTL) or you manually executed the kill command. | No further action is required. If the task was terminated due to TTL, consider adjusting the allowed duration if needed.
ZOMBIE | The task has stopped unexpectedly because the computational resources were shut down. | The task was still in the SUBMITTED status when the machine group was shut down. Since the task hadn’t started, there are no outputs or logs. | Restart the task by submitting it to a new machine group or ensure the machine group remains active for future tasks.
SPOT PREEMPTED | The spot instances running your task were terminated, and the task has been paused. | Spot instances were reclaimed by the cloud provider. If auto-resubmission is enabled, the task is re-queued and goes back to SUBMITTED status when new resources become available. Otherwise, it stays in this state. | If cost savings are a priority, spot instances are a good choice. For uninterrupted tasks, consider switching to on-demand instances. If the task is stuck, enable auto-resubmission or manually resubmit it.
FAILED | The task failed due to an error in the simulator, likely caused by incorrect input configurations or an internal simulator issue. | The simulator command returned a non-zero status code, indicating failure. Error details are available only if the user chose task.wait. | Inspect the logs (stderr and stdout) for more information. If you didn’t use task.wait, resubmit the task with better input configurations.
COMPUTATION ENDED | The simulation's commands have ended. | The simulation commands have ran. Results will now be uploaded to the user's storage. | The task's results will be uploaded to the user's storage and the task will go to a final status.
EXECUTOR TERMINATED | The executor running your task was terminated due to internal reasons from the provider. | The cloud provider terminated the executor unexpectedly. If auto-resubmission is enabled, the task is re-queued and returns to SUBMITTED status. Otherwise, it remains in this state. | Enable auto-resubmission to avoid delays, or manually resubmit the task. For more stability, consider using on-demand instances if spot instances were in use.
EXECUTOR TERMINATED BY USER | The task’s executor was terminated by you while it was running. | The machine group was terminated after the task had started, but no outputs or logs are available because they are only saved once the task completes. | If outputs are needed, avoid terminating the machine group mid-task. Resubmit the task and let it finish to retrieve logs and results.
EXECUTOR TERMINATED TTL EXCEEDED | The task’s executor was terminated because it exceeded its time-to-live (TTL) limit. | The Machine Group’s TTL, defined by your quotas, was reached. This helps manage costs by ensuring resources don’t run longer than expected. | Start a new Machine Group and resubmit the task. Alternatively, increase your quotas or manually set a TTL if necessary to avoid future interruptions.
EXECUTOR FAILED | The task failed due to an executor error, such as low disk space. | An exception in the API, most commonly due to insufficient disk space, caused the task to fail. Error details are only shown if you used task.wait. | Check your machine group’s disk space or other resource limits. Resubmit the task after resolving the issue, and consider using task.wait for future tasks to see real-time error messages.
TTL EXCEEDED | The task exceeded its configured time-to-live (TTL) and was automatically stopped. | The task ran longer than the TTL defined in your quotas, which helps control costs by limiting how long a task can use shared resources. You can also set a TTL manually for better control. | Resubmit the task and consider adjusting the TTL if more computation time is needed. Alternatively, increase your quotas to avoid future interruptions.
