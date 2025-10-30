# task-runner

## inductiva task-runner - CLI interface

```default
inductiva task-runner [-h] {launch,remove} ...
```

Task-Runner management utilities.

A TaskRunner is a local process that manages and executes Inductiva tasks. With the Inductiva CLI, you can launch runners locally, so you can use your own resources to run your tasks. This allows you to switch between local and remote execution, and balance between cost (zero compute cost if you are running Tasks on your local machines) and performance (access to high-performance machines on the cloud via Inductiva). 

The `inductiva task-runner` command allows you to manage your local task-runners on the Inductiva platform. It provides utilities including launching and terminating task-runners.

### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

---

### inductiva task-runner launch

```default
inductiva task-runner launch [-h] [-w [WATCH]] [--provider {local,gcp}] [--hostname HOSTNAME]
                             [--detach] [--zone ZONE] [--machine-type MACHINE_TYPE] [--spot]
                             [--max-idle-time MAX_IDLE_TIME] [--disk-size DISK_SIZE]
                             machine_group_name
```

The `inductiva task-runner launch` command provides a way to launch a Task-Runner on different providers.

Providers:
  local  - Launch locally using Docker (default)
  gcp    - Launch on Google Cloud Platform

#### Positional Arguments

* [**`machine_group_name`**]() - Name of the machine group to launch the Task-Runner. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-w`**]() `WATCH`, [**`--watch`**]() `WATCH` - Prompt the command every N seconds. (default: `None`)
* [**`--provider`**]() `PROVIDER`, [**`-p`**]() `PROVIDER` - Provider to use for launching the task-runner (default: local). (default: `local`)
* [**`--hostname`**]() `HOSTNAME`, [**`-ho`**]() `HOSTNAME` - Hostname of the Task-Runner. (default: `None`)
* [**`--detach`**](), [**`-d`**]() - Run the task-runner in the background (local only).

---

#### Options

* [**`--zone`**]() `ZONE`, [**`-z`**]() `ZONE` - GCP zone where the VM will be created (default: europe-west1-b). (default: `europe-west1-b`)
* [**`--machine-type`**]() `MACHINE_TYPE`, [**`-t`**]() `MACHINE_TYPE` - GCP machine type (default: c2d-standard-8). (default: `c2d-standard-8`)
* [**`--spot`**]() - Use spot instance.
* [**`--max-idle-time`**]() `MAX_IDLE_TIME` - Idle time in minutes before auto-termination (default: 3). (default: `None`)
* [**`--disk-size`**]() `DISK_SIZE`, [**`-ds`**]() `DISK_SIZE` - Disk size in GB (default: 10). (default: `10`)

---

### inductiva task-runner remove

```default
inductiva task-runner remove [-h] [-w [WATCH]]
```

The `inductiva task-runner remove` command terminates and removes a running task-runner.

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-w`**]() `WATCH`, [**`--watch`**]() `WATCH` - Prompt the command every N seconds. (default: `None`)

---
::docsbannersmall
::
