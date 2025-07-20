# Tasks

The Inductiva API provides complete task management for tracking, monitoring, and controlling your simulation workloads. Tasks represent the complete lifecycle of your simulations, from submission to completion.


## Get started
Learn about _Tasks_ — the execution unit in the Inductiva API that represents your simulation workloads. Tasks provide real-time monitoring, status tracking, and result management for all your computational work.

| **[Tasks →](tasks.md)** | **[Task Execution →](#task-execution)** | **[Task Lifecycle →](#task-lifecycle)** |
|---|---|---|
| Learn about Tasks and how they are automatically generated when you submit simulations | Learn about task execution flow with synchronous and asynchronous approaches | Understand the task's status state machine, from submission through completion |

## Features
✓ **Automatic task creation** Every simulation submission creates a unique Task with comprehensive tracking. No manual setup required — just run your simulation and get instant task management.

✓ **Asynchronous execution** Submit multiple simulations without blocking your workflow. Tasks run in parallel while you continue developing, with full control over when to wait for completion.

✓ **Real-time monitoring** Track task progress through detailed status updates and lifecycle stages. Monitor resource usage, execution time, and identify bottlenecks in real-time.

✓ **Persistent task metadata** Store simulation parameters and custom data directly with your tasks. Eliminate external tracking systems by keeping all relevant information integrated with your computational workloads.

✓ **Robust error handling** Automatic recovery from infrastructure interruptions including spot instance reclamation and machine termination. Failed tasks provide detailed diagnostics for quick troubleshooting.

✓ **Session independence** Recreate task objects across different sessions using task IDs. Continue monitoring and managing tasks even after restarting your development environment.