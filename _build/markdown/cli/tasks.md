# tasks

## inductiva tasks - CLI interface

```default
inductiva tasks [-h] {download,info,kill,last-modified-file,list,list-files,tail,top} ...
```

Task management utilities.

The `inductiva tasks` command provides tools to manage your tasks on  the Inductiva platform. It includes utilities for monitoring task  progress, viewing details, and terminating tasks when needed.

### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

---

### inductiva tasks download

```default
inductiva tasks download [-h] [--filenames [FILENAMES ...]] [--dir DIR] [--input] [--output]
                         [task_id ...]
```

Download the input and/or output files of tasks with the given ID(s).

By default, all output files of the provided task are downloaded to a local directory at `inductiva_output/<TASK_ID>/outputs`, relative to your current working directory.

To download specific files, use the `--filenames` option and provide a list of filenames.

You can specify multiple task IDs to download files from several tasks at once.                                    

The target directory can be set with the `--dir` option. Files will be saved into a subdirectory named after each task ID.

Use the `--output` (`-o`) option to download output files, and the
`--input` (`-i`) option to download input files. If neither is specified, only the output files are downloaded by default.

#### Positional Arguments

* [**`task_id`**]() - ID(s) of the task(s) to download. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`--filenames`**]() `FILENAMES` - Names of the files to download. (default: `None`)
* [**`--dir`**]() `DIR` - Path to the directory where input/output files will be downloaded. (default: `None`)
* [**`--input`**](), [**`-i`**]() - Option to download input files.
* [**`--output`**](), [**`-o`**]() - Option to download output files.

---

### inductiva tasks info

```default
inductiva tasks info [-h] [-w [WATCH]] id
```

The `inductiva tasks info` command provides an extensive list of information about a task.

#### Positional Arguments

* [**`id`**]() - ID of the task to get information about. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-w`**]() `WATCH`, [**`--watch`**]() `WATCH` - Prompt the command every N seconds. (default: `None`)

---

### inductiva tasks kill

```default
inductiva tasks kill [-h] [-w WAIT_TIMEOUT] [-y] [--all] [id ...]
```

The `inductiva tasks kill` command terminates the specified tasks on the platform. You can terminate multiple tasks by passing multiple task IDs. To confirm termination without a prompt, use the `-y` or `--yes` option. If you provide `-w` or `--wait-timeout`, the system does not confirm whether the termination was successful.

Note: The `inductiva tasks kill` command does not stop the machine where the task is running. It only terminates the task itself, leaving the computational resources active and available to run other tasks.

#### Positional Arguments

* [**`id`**]() - ID(s) of the task(s) to kill. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-w`**]() `WAIT_TIMEOUT`, [**`--wait-timeout`**]() `WAIT_TIMEOUT` - Number of seconds to wait for the kill command. If not provided, the system sends the request without waiting a response. (default: `None`)
* [**`-y`**](), [**`--yes`**]() - Skip kill confirmation.
* [**`--all`**]() - Kill all running tasks.

#### Examples

```bash

$ inductiva tasks kill cmvsc9qhz5iy86f6pef8uyxqt
You are about to kill the following tasks:
    - cmvsc9qhz5iy86f6pef8uyxqt 
Are you sure you want to proceed (y/[N])? y
Successfully sent kill request for task cmvsc9qhz5iy86f6pef8uyxqt.
```

---

### inductiva tasks last-modified-file

```default
inductiva tasks last-modified-file [-h] id
```

The `inductiva tasks last-modified-file` command displays the last modified file for the simulation. It also shows the last time the file was modified and how long it has been since that.

#### Positional Arguments

* [**`id`**]() - The ID of the task to get the last modifed file. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

#### Examples

```bash

$ inductiva tasks last-modified-file qpusar8bch509k56g1hvv5yxk

Most Recent File: /mnt/disks/task-runner-data/workdir/qpusar8bch509k56g1hvv5yxk/output/artifacts/stdout.txt
Modification Time: 2025-04-03 12:58:49
Current Time on Machine: 2025-04-03 12:58:50

Time Since Last Modification: 0:00:01
```

---

### inductiva tasks list

```default
inductiva tasks list [-h] [-w [WATCH]] [-n LAST_N | -a | -i TASK_ID [TASK_ID ...]]
                     [-p PROJECT_NAME]
```

The `inductiva tasks list` command provides an overview of your tasks on the platform. It lists the most recent tasks or a specific task by ID. You can control the number of tasks listed with the `-n` or `--last-n` options.

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-w`**]() `WATCH`, [**`--watch`**]() `WATCH` - Prompt the command every N seconds. (default: `None`)
* [**`-n`**]() `LAST_N`, [**`--last-n`**]() `LAST_N` - List most recent tasks. (default: `10`)
* [**`-a`**](), [**`--all`**]() - List all tasks.
* [**`-i`**]() `TASK_ID`, [**`--id`**]() `TASK_ID` - List a task with a specific ID. (default: `None`)
* [**`-p`**]() `PROJECT_NAME`, [**`--project-name`**]() `PROJECT_NAME` - List the tasks of a project. (default: `None`)

#### Examples

```bash

# List the last 4 tasks with their details:
$ inductiva tasks list -n 4
ID                              SIMULATOR          STATUS         SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
jxwt0rm8s8xspdfcegtgkkana       splishsplash       started        08 Feb, 13:25:49       08 Feb, 13:26:04       *0:00:05                 c2-standard-4
n0zcac8rmw7xydbis3m407kb4       splishsplash       started        08 Feb, 13:25:48       08 Feb, 13:26:03       *0:00:07                 c2-standard-4
8nmpn4h99nyfpo4da9jw2405q       splishsplash       started        08 Feb, 13:25:47       08 Feb, 13:26:02       *0:00:09                 c2-standard-4
so6i93pi74b89rndircubp3v2       splishsplash       started        08 Feb, 13:25:47       08 Feb, 13:26:02       *0:00:10                 c2-standard-4
```

---

### inductiva tasks list-files

```default
inductiva tasks list-files [-h] [-w [WATCH]] id
```

The `inductiva tasks list-files` command lists the contents of a task’s working directory while the task is in progress. (Experimental)

#### Positional Arguments

* [**`id`**]() - ID of the task to list files. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-w`**]() `WATCH`, [**`--watch`**]() `WATCH` - Prompt the command every N seconds. (default: `None`)

---

### inductiva tasks tail

```default
inductiva tasks tail [-h] [--lines LINES] [--follow] id filename [filename ...]
```

The `inductiva tasks tail` shows the last lines of a file in the directory of a task that is running. (Experimental)

#### Positional Arguments

* [**`id`**]() - ID of the task to list directories. (default: `None`)
* [**`filename`**]() - File to tail. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`--lines`**]() `LINES`, [**`-l`**]() `LINES` - Number of lines to show. (default: `10`)
* [**`--follow`**](), [**`-f`**]() - Keep the file open and show new lines.

---

### inductiva tasks top

```default
inductiva tasks top [-h] id
```

The `inductiva tasks top` command streams the output of the `top -b -H -n 1` command executed on the machine where the task is running. This allows real-time monitoring of system processes and resource usage for the task.

#### Positional Arguments

* [**`id`**]() - The ID of the task to monitor. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

#### Examples

```bash

$ inductiva tasks top qpusar8bch509k56g1hvv5yxk
top - 12:00:15 up 18 min,  0 users,  load average: 1.14, 0.99, 0.58
Threads: 226 total,   2 running, 224 sleeping,   0 stopped,   0 zombie
%Cpu(s): 24.2 us,  1.5 sy,  0.0 ni, 72.7 id,  0.0 wa,  0.0 hi,  1.5 si,  0.0 st
MiB Mem :  16008.2 total,  12976.4 free,   1057.1 used,   1974.7 buff/cache
MiB Swap:      0.0 total,      0.0 free,      0.0 used.  14656.3 avail Mem 

    PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND
    1469 task-ru+  20   0  894208 711108  36128 R  99.9   4.3   9:56.46 d_hydro+
    1557 task-ru+  20   0    9016   3812   3140 R   6.2   0.0   0:00.01 top
        1 root      20   0  165128  10828   7912 S   0.0   0.1   0:01.22 systemd
        2 root      20   0       0      0      0 S   0.0   0.0   0:00.00 kthreadd
        3 root       0 -20       0      0      0 I   0.0   0.0   0:00.00 rcu_gp
        4 root       0 -20       0      0      0 I   0.0   0.0   0:00.00 rcu_par+
...
```

---
::docsbannersmall
::
