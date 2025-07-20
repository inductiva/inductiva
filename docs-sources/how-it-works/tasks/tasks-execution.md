# Task Execution & Management

Tasks run **asynchronously by default**, allowing you to submit simulations and continue working without blocking your code. This design enables efficient batch processing, parallel execution across multiple resources, and resilient handling of long-running computations.

## Core Concepts

### Asynchronous Execution
When you submit a task by doing `simulator.run()`, it immediately returns a `Task` object while the simulation runs in the background. This allows you to:

- Submit multiple simulations in parallel
- Continue local work while computations run remotely
- Handle task completion at your own pace
- Maintain fault tolerance for long-running jobs

### Task Lifecycle
Each task progresses through several main states:

| Status | Description | Next States |
|--------|-------------|-------------|
| `pending-input` | Task is waiting for input files | `in-queue` |
| `in-queue` | Task is queued in machine waiting for it to become available | `submitted` |
| `submitted` | Task is queued and waiting for resources | `started`, `failed` |
| `started` | Task is actively running on compute resources | `success`, `failed` |
| `success` | Task finished successfully | - |
| `failed` | Task encountered an error | - |
| `killed` | Task was manually terminated | - |

````{eval-rst}
.. seealso::
    For complete Task Lifecycle documentation, see our `Task Lifecycle <https://inductiva.ai/guides/how-it-works/intro/tasks-lifecycle>`_ guide
````

## Execution Control

### Blocking Execution with `wait()`

Convert asynchronous tasks into blocking operations by using the `wait` method. This method will block the call, forcing your script to wait until the task reaches a terminal status to continue its execution.

This is particularly useful when you want to wait for the completion of a task before proceeding with some other operation (e.g., downloading the output files):

```python
# Submit task (non-blocking)
task = simulator.run(input_dir="config/", on=machine_group)

# Block the call until the simulation completes
task.wait()  # <- The remote simulation WILL NOT DIE if the local session is
             #    interrupted while waiting for the wait() call to return

# Now safe to download results
task.download_outputs()
```

> Note: The `wait()` method is resilient to local interruptions. Your remote simulation will continue running even if your Python script is terminated.

### Parallel Execution Workflow

Leverage asynchronous execution to run multiple simulations simultaneously:

```python
import inductiva

# Initialize resources
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()
simulator = inductiva.simulators.SplishSplash()

# Submit batch of tasks
tasks = []
configs = ["wind_10ms", "wind_15ms", "wind_20ms", "wind_25ms"]

for config in configs:
    task = simulator.run(
        input_dir=f"scenarios/{config}",
        on=machine_group
    )
    tasks.append(task)
    print(f"✓ Submitted {config}: {task.id}")

print(f"Submitted {len(tasks)} tasks in parallel")

# Process results as they complete
for i, task in enumerate(tasks):
    print(f"⏳ Waiting for task {i+1}/{len(tasks)}...")
    task.wait()
    
    if task.status == "completed":
        print(f"✅ {configs[i]} completed successfully")
        task.download_outputs(f"results/{configs[i]}")
    else:
        print(f"❌ {configs[i]} failed: {task.status}")

machine_group.terminate()
```

## Task Management

The [Python Client](https://inductiva.ai/guides/api-functions/api/inductiva.tasks) provides several methods for managing tasks.

| Method | Description | Example |
|--------|-------------|---------|
| `inductiva.tasks.get()` | Retrieve recent tasks with optional filtering | `inductiva.tasks.get(last_n=10, status="completed")` |
| `inductiva.tasks.get_all()` | Fetch all tasks (with pagination) | `inductiva.tasks.get_all(project="fluid-dynamics")` |
| `inductiva.tasks.get_tasks()` | Flexible task retrieval with smart defaults | `inductiva.tasks.get_tasks(last_n=5, status="started")` |
| `task.wait()` | Block until task reaches terminal state | `task.wait()` |
| `task.kill()` | Terminate a running or queued task | `task.kill()` |
| `task.get_info()` | Fetch detailed task metadata | `task.get_info()` |
| `task.get_position_in_queue()` | Check queue position for pending tasks | `task.get_position_in_queue()` |
| `task.download_outputs()` | Download simulation results | `task.download_outputs("./results")` |

Inductiva also provides [CLI methods](https://inductiva.ai/guides/api-functions/cli/tasks) for managing your tasks.

| Command | Description | Example |
| :--- | :--- | :--- |
| `inductiva tasks list` | List all tasks with status and details | `inductiva tasks list -n 4` |
| `inductiva tasks info` | Get detailed information about a task | `inductiva tasks info <task-id>` |
| `inductiva tasks download` | Download task outputs to a local machine | `inductiva tasks download <task-id>` |
| `inductiva tasks kill` | Terminate a running or queued task | `inductiva tasks kill <task-id>` |
| `inductiva tasks list-files` | Show current files of a running task | `inductiva tasks list-files --id <task-id>` |
| `inductiva tasks tail` | Display last lines of task output | `inductiva tasks tail --id <task-id> -L 50` |
| `inductiva tasks top` | Show system processes on the task machine | `inductiva tasks top <task-id>` |
| `inductiva tasks last-modified-file` | Show the most recently modified file | `inductiva tasks last-modified-file <task-id>` |
| `inductiva logs` | Stream task logs in real-time | `inductiva logs <task-id>` |