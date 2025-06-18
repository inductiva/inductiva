# Command Line Interface (CLI)

Welcome to the **Inductiva CLI**, the command-line tool that enables
users to interact with the Inductiva platform directly from the terminal window.

## Available Commands

The CLI provides the following commands:

### General Commands

- **`auth`** - Manage authentication and API tokens.
- **`user`** - View and update user settings.

### Task Management

- **`tasks`** - List, start, stop, and monitor tasks.
- **`task-runner`** - Manage the execution of local Task Runners.
- **`logs`** - View logs of running tasks.

### Project and Storage Management

- **`projects`** - Manage projects.
- **`storage`** - Upload and manage remote storage files.

### Computational Resources

- **`resources`** - View and allocate computational resources.
- **`simulators`** - List available simulators.

## Set Up & Authentication

The CLI is automatically installed when you install
Inductiva's Python Client. However, before using the CLI,
you need to authenticate with your Inductiva account.
From your terminal window type:

```sh
inductiva auth login
```

When prompted, copy paste your personal API key that is availalbe
from your [User Account] (<https://console.inductiva.ai/account/profile>)
page in the Web Console.

## Need Help?

Use the `--help` or `-h` flag with **any command** for more details:

```sh
inductiva tasks --help
```
