# task-runner

The `inductiva task-runner` command allows you to manage local
Task-Runners on the Inductiva platform. Task-runners are responsible for executing 
tasks **locally** instead of on remote machines.

## Usage

```sh
inductiva task-runner [-h] {launch,remove} ...
```

### Description
A TaskRunner is a local process that manages and executes Inductiva Tasks.
With the Inductiva CLI, you can launch runners locally, so you can use your
own resources to run your Tasks. This allows you to switch between local and
remote execution, and balance between cost (zero compute cost if you are 
running Tasks on your local machines) and performance (access to high-performance
machines on the cloud via Inductiva). 

## Options

- **`-h, --help`** → Show help message and exit.

## Available Subcommands

### `launch`
Launch a new task-runner instance.

```sh
inductiva task-runner launch [-h] [-w [WATCH]] [--hostname HOSTNAME] [--detach] machine_group_name
```

#### Options for `launch`:
- **`-h, --help`** → Show help message and exit.
- **`-w [WATCH]`** → Enable watch mode to monitor the task-runner's status.
- **`--hostname HOSTNAME`** → Specify the hostname for the task-runner.
- **`--detach`** → Run the task-runner in detached mode.
- **`machine_group_name`** → The name of the machine group where the task-runner will be launched.

### `remove`
Terminate and remove an existing task-runner.

```sh
inductiva task-runner remove
```

## Example Usage

### Start a task-runner:
```sh
inductiva task-runner launch my_machine_group
```

### Start a detached task-runner:
```sh
inductiva task-runner launch my_machine_group --detach
```

### Remove a task-runner:
```sh
inductiva task-runner remove
```

## Need Help?
Run the following command for more details:

```sh
inductiva task-runner --help
```
