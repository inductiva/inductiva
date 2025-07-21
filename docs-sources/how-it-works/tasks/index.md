# Tasks

Once you've launched a simulation how do you track it? The answer is **Tasks**. Inductiva API creates a unique Task for every simulation run, providing a dedicated object to monitor, manage, and retrieve results related to the computation from start to finish.


## Get started
Learn about _Tasks_ — the execution unit in the Inductiva API that represents your simulation workloads. Tasks provide real-time monitoring, status tracking, and result management for your computational work.

| **[Tasks →](tasks.md)** | **[Task Execution →](tasks-execution.md)** | **[Task Lifecycle →](tasks-lifecycle.md)** |
|---|---|---|
| Learn about Tasks and how they are automatically generated when you submit simulations | Control a task execution flow with synchronous and asynchronous approaches | Understand the task's status state machine, from submission through completion |

## Features
✓ **Automatic task creation** Every simulation submission creates a unique Task for managing that simulation. No manual setup required — just run your simulation and get instant task management.

✓ **Asynchronous execution** Submit multiple simulations without blocking your workflow. Tasks run in parallel while you continue developing, with full control over when to wait for completion.

✓ **Real-time monitoring** Track task progress through status update. Monitor resource usage, execution time, and identify system bottlenecks in real-time.

✓ **Persistent task metadata** Store simulation parameters and custom data directly with your tasks. Eliminate external tracking systems by keeping all relevant information integrated.

✓ **Session independence** Recreate task objects across different sessions using task IDs. Continue monitoring and managing tasks even after restarting your development environment.