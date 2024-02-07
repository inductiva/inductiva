# Command Line Interface (CLI)

Walking side by side with the Python client, the Inductiva CLI provides 
a quick and simple interface to interact and perform the most common tasks from
the terminal. The CLI allows users to quickly manage their tasks, resources and/or
storage, while their simulations are still running remotely.

## Overview of the CLI

The CLI beauty is that its functionalities are a terminal command distance away.
By doing,
```bash
inductiva
```

you will obtain an overview of the available commands and their descriptions.

Here, we highlight and go deeper into the four main functionalities of the CLI:
- [`resources`](./resources) give management capabilities over the resources;
- [`tasks`](./tasks#task-management) empower the quick management of the simulations;
- [`logs`](./tasks#streaming-logs-of-a-task) allow users to stream the logs of their running tasks;
- [`storage`](./storage) provides an overview of the user's remote storage.

Whenever you need help on a specific command, you can always use the `--help` or
`-h` flag to further understand its usage.

```{toctree}
---
hidden: true
---
resources
tasks
storage
```