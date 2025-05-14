# tasks
The Inductiva API allows you to run multiple simulations remotely, and for every
simulation you submit to the API, you generate a separate Task. 
The `inductiva tasks` command allows you to inspect and manage all the tasks
you are run on Inductiva.


## Usage

```sh
inductiva tasks [-h] {download,info,kill,list,list-files,tail} ...
```

### Options
- **`-h, --help`** → Show help message and exit.

## Available Subcommands

### `download`
Downloads the output of completed tasks to your local machine. 

```sh
inductiva tasks download TASK_ID
```
By default, the outputs of  the task (multiple files) will be
stored in a local directory named `output_dir/tasks_id` inside
your project directory. You can specify an alternative the directory 
for storing the files with `--output_dir`. 

You can also download the outputs from multiple tasks:
```bash
$ inductiva tasks download task_1_id task_2_id
```

You can download specific files from the output of the task by 
setting the `--filenames` flag.


### `info`
Retrieve detailed information about a specific task.

```sh
inductiva tasks info TASK_ID
```

### `kill`
Kill a running task. 
```sh
inductiva tasks kill TASK_ID
```
The kill subcommand will prompt you for a confirmation so that you do not
accidentally kill an important task. Note, however, that the kill subcommand
does not stop the machine where the task was running: it just terminates the
task, leaving the computational resources where it was running ready for
taking other tasks.

### `last-modified-file`

Displays the last modified file of a task.

```sh
inductiva tasks last-modified-file TASK_ID
```

### `list`
List all tasks associated with your account.

```sh
inductiva tasks list
```

### `list-files`
Show the current files of a running task.

```sh
inductiva tasks list-files --id TASK_ID
```

### `tail`
Display the last lines of a task’s file output.

```sh
inductiva tasks tail --id TASK_ID -l NUMBER_OF_LINES
```

### `top`

Displays the output of the `top` command from the machine running the task.

```sh
inductiva tasks top TASK_ID
```

## Example Usage

List all tasks:

```sh
inductiva tasks list
```

List the last four tasks with their details

```sh
$ inductiva tasks list -n 4
       ID                              SIMULATOR          STATUS         SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
       jxwt0rm8s8xspdfcegtgkkana       splishsplash       started        08 Feb, 13:25:49       08 Feb, 13:26:04       *0:00:05                 c2-standard-4
       n0zcac8rmw7xydbis3m407kb4       splishsplash       started        08 Feb, 13:25:48       08 Feb, 13:26:03       *0:00:07                 c2-standard-4
       8nmpn4h99nyfpo4da9jw2405q       splishsplash       started        08 Feb, 13:25:47       08 Feb, 13:26:02       *0:00:09                 c2-standard-4
       so6i93pi74b89rndircubp3v2       splishsplash       started        08 Feb, 13:25:47       08 Feb, 13:26:02       *0:00:10                 c2-standard-4
```

```sh
inductiva tasks info --id TASK_ID
```

Terminate a task:

```sh
inductiva tasks kill --id TASK_ID
```
For a specific task, we will get something like:
```sh
$ inductiva tasks kill cmvsc9qhz5iy86f6pef8uyxqt
You are about to kill the following tasks:
  - cmvsc9qhz5iy86f6pef8uyxqt 
Are you sure you want to proceed (y/[N])? y
Successfully sent kill request for task cmvsc9qhz5iy86f6pef8uyxqt.
```

Get the last modified file of a task:

```sh
inductiva tasks last-modified-file qpusar8bch509k56g1hvv5yxk

Most Recent File: /mnt/disks/task-runner-data/workdir/qpusar8bch509k56g1hvv5yxk/output/artifacts/stdout.txt
Modification Time: 2025-04-03 12:58:49
Current Time on Machine: 2025-04-03 12:58:50

Time Since Last Modification: 0:00:01
```

Display the processes running on the machine where the task is running:

```sh
inductiva tasks top qpusar8bch509k56g1hvv5yxk
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

## Need Help?
Run the following command for more details:

```sh
inductiva tasks --help
```