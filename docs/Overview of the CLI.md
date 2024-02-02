## Command Line Interface (CLI)

Walking side by side with the Python client, the Inductiva CLI provides 
a quick and simple interface to interact and perform the most common tasks from
the terminal. The CLI allows users to quickly manage their tasks, resources and/or
storage, while their simulations are still running remotely.

### Overview of the CLI

The CLI beauty is that its functionalities are a terminal command distance away.
By doing,
```bash
inductiva
```

you will obtain an overview of the available commands and their descriptions.

Here, we highlight and go deeper on the four main functionalities of the CLI:
- [`resources`](#resource-management-in-the-cli) gives management capabilities over the resources;
- [`tasks`](#task-management) empowers the quick management of the simulations;
- [`logs`](#streaming-logs-of-a-task) allows users to stream the logs of their running tasks;
- [`storage`](#remote-storage-overview) provides an overview of the user's remote storage.


Whenever you need help on a specific command, you can always use the `--help` or
`-h` flag to further understand its usage.

#### Resource Management in the CLI

As you might know by now, Inductiva API provides a simple way to [launch
dedicated resources]() where you can run your simulations. Before launching any
resources, users can use the CLI to gather more information about the right
resources to launch with the `available` and `cost` subcommands.

With them, user can list all the available machine types together with details,
```bash
$ inductiva resources available
Available machine types

machine-type: [cores-available]
c2-standard-: [4, 8, 16, 30, 60]
c3-standard-: [4, 8, 22, 44, 88, 176]
c2d-standard-: [2, 4, 8, 16, 32, 56, 112]
c2d-highcpu-: [2, 4, 8, 16, 32, 56, 112]
e2-standard-: [2, 4, 8, 16, 32]
n2-standard-: [2, 4, 8, 16, 32, 48, 64, 80, 96, 128]
n2d-standard-: [2, 4, 8, 16, 32, 48, 64, 80, 96, 128, 224]
n1-standard-: [1, 2, 4, 8, 16, 32, 64, 96]
e2-highcpu-: [2, 4, 8, 16, 32]

 E.g. of machine-type: c2-standard-8
```

and thereafter, use one of the machine types to get an estimate of the cost
per hour to use it:
```bash
$ inductiva resources cost c2-standard-8 --spot -n 4
Estimated total cost (per machine): 0.445 (0.111) $/h.
```

When you have already decided and launch a few computational resources, you can
use the list subcommand to get an overview of the resources you have running:
```bash
$ inductiva resources list
Active Resources:
Name                           Machine Type    Elastic    Type        # machines    Disk Size in GB  Spot    Started at (UTC)
-----------------------------  --------------  ---------  --------  ------------  -----------------  ------  ------------------
api-abs7c8hwccu0v6cppo9nc40z6  c2-standard-30  False      standard             5                 70  True    01 Feb, 23:06:52
api-mym74u9i9xppi9ngzze3xueae  c2-standard-8   False           mpi             4                 70  True    02 Feb, 00:29:53
```

Finally, the CLI also allows to terminate the resources that are no longer required with
the `terminate` subcommand. Users can either choose a specific resource with the
`--name` flag or terminate all the resources with the `--all` flag. Any of the steps
require confirmation from the user before proceeding.
```bash
$ inductiva resources terminate --all
Confirm the termination of all active machines? (y/[N]) y
Terminating all active computational resources.
Terminating MachineGroup(name="api-abs7c8hwccu0v6cppo9nc40z6"). This may take a few minutes.
Machine Group api-abs7c8hwccu0v6cppo9nc40z6 with c2-standard-30 machines successfully terminated in 0:01:07.
Terminating MachineGroup(name="api-mym74u9i9xppi9ngzze3xueae"). This may take a few minutes.
Machine Group api-mym74u9i9xppi9ngzze3xueae with c2-standard-8 machines successfully terminated in 0:01:12.
```

### Task Management

The power of Inductiva API lies on the ability to launch multiple simulations remotely.
As an easy to track them, the CLI provides a `list` subcommand that provides information
about the last tasks launched by the user. If you wish, you can also just track
a specific task via its `--id`. Examples of usage:

```bash
$ inductiva tasks list -n 4
ID                         Simulator    Status    Submitted         Started           Computation Time    Resource Type
-------------------------  -----------  --------  ----------------  ----------------  ------------------  ---------------
s9evg758ieteqzvg7p6hnuqcw  swash        started   02 Feb, 00:52:10  02 Feb, 00:52:14  *0:00:19            c2-standard-4
tc7cwuer45kfzuw8t93r6dxa8  swash        success   01 Feb, 23:52:44  01 Feb, 23:52:49  0:25:26             c2-standard-4
hzgk5ngzk28a39qa7mesv0snk  swash        success   01 Feb, 23:45:19  01 Feb, 23:45:19  0:25:51             c2-standard-4
mjnb8c7i8bfppgmu2y1zd1o7f  swash        success   01 Feb, 23:24:19  01 Feb, 23:25:08  0:09:43             c2-standard-30
```

One cool tip is to use a [watch method](https://www.geeksforgeeks.org/watch-command-in-linux-with-examples/) to recurrently monitor the task information.
For example in Linux or Mac, one can monitor the task information every 10 seconds with:
```bash
$ watch -n 10 inductiva tasks -id tc7cwuer45kfzuw8t93r6dxa8
Every 10.0s: inductiva tasks list -id s9evg758ieteqzvg7p6hnuqcw                                      

ID                         Simulator    Status    Submitted         Started           Computation Time    Resource Type
-------------------------  -----------  --------  ----------------  ----------------  ------------------  ---------------
s9evg758ieteqzvg7p6hnuqcw  swash        started   02 Feb, 00:52:10  02 Feb, 00:52:14  *0:01:27            c2-standard-4
```

In case, you are not satisfied with the results of a running simulation, you can
always `kill` it and free up the computational resources for other tasks. E.g.:

```bash
$ inductiva tasks kill -id s9evg758ieteqzvg7p6hnuqcw
Do you confirm you want to kill 1 tasks (y/[N])?y
Successfully sent kill request for task s9evg758ieteqzvg7p6hnuqcw.
> Consider tracking the task status by calling task.get_status()
> or via the CLI using 'inductiva tasks -id s9evg758ieteqzvg7p6hnuqcw'.
```

### Streaming Logs of a Task

While a task is running, one of the main concerns for scientists and engineers is to understand if the simulation is progressing as expected and/or if there are any errors.

To help with this, the CLI provides a `logs` subcommand that allows users to stream the logs of a running task in real time. In this way, running a simulation in a powerful remote machine will
feel exactly like running in your local computer. Example:

```bash
$ inductiva logs f0bnqgf4fcr4asgi4e21tcsqa
Websocket connected
Opening socket connection to logs of task f0bnqgf4fcr4asgi4e21tcsqa ...
timestep: 0.0226
t/T: 12.2
breaking: 0
Fi_iter: 4 Final_residual: 3.73e-09  Fi_time: 0.00174
breaking: 0
Fi_iter: 4 Final_residual: 3.08e-09  Fi_time: 0.00172
breaking: 0
Fi_iter: 5 Final_residual: 1.68e-11  Fi_time: 0.00167
umax: 0.102
vmax: 0
wmax: 0.0941
dt: 0.0226
wavegentime: 2.49e-05
reinitime: 0
gctime: 0.000539         average gctime: 0.0005
Xtime: 0.000126  average Xtime: 0.000641
total time: 7.62451   average time: 0.00878
timer per step: 0.00885
------------------------------------
869
simtime: 19.6
timestep: 0.0226
t/T: 12.2
breaking: 0
Fi_iter: 4 Final_residual: 2.68e-09  Fi_time: 0.00187
breaking: 0
Fi_iter: 4 Final_residual: 5.84e-09  Fi_time: 0.00187
breaking: 0
Fi_iter: 4 Final_residual: 9.24e-09  Fi_time: 0.00145
umax: 0.101
vmax: 0
wmax: 0.0984
dt: 0.0226
wavegentime: 2.06e-05
reinitime: 0
gctime: 0.000429         average gctime: 0.0005
Xtime: 0.000979  average Xtime: 0.000641
total time: 7.63317   average time: 0.00878
timer per step: 0.00865
------------------------------------
870
simtime: 19.6
timestep: 0.0226
t/T: 12.2
breaking: 0
Fi_iter: 4 Final_residual: 2.3e-09  Fi_time: 0.00165
breaking: 0
Fi_iter: 5 Final_residual: 5.82e-12  Fi_time: 0.00194
breaking: 0
Fi_iter: 4 Final_residual: 3.39e-09  Fi_time: 0.00143
umax: 0.101
vmax: 0
wmax: 0.109
dt: 0.0226
wavegentime: 2.39e-05
reinitime: 0
gctime: 0.00056  average gctime: 0.0005
Xtime: 0.00019   average Xtime: 0.000641
total time: 7.64189   average time: 0.00878
timer per step: 0.00873
------------------------------------
```

This combined with the kill method, allows a quick iteration to stabilize the simulation and achieve the desired configuration.

### Remote Storage overview

Finally, when the simulation finishes the results are saved in the user's remote bucket. 

Hence, the CLI allows users to connect to their remote bucket where all the data of simulations live and manage it as they wish.

To start, they can explore the storage use at any time with:

```bash
$ inductiva storage size
Total user's remote storage in use: 0.46 GB
```

Thereafter, the contents of the storage can be listed as follows:
```bash
$ inductiva storage list --max-results 10 --order-by size --sort-order desc
Name                        Size      Creation Time
--------------------------  --------  ----------------
hodbisrxjxhdbkknv60xmy6ti/  87.55 MB  01 Feb, 14:55:43
m4d487kf46n4pp868qddzn67m/  54.48 MB  01 Feb, 16:49:36
ibeg639dv96yo2kx18wbtxsfi/  51.93 MB  01 Feb, 14:57:25
2oqv8tyfa1dubeq63z4g5e2zn/  36.67 MB  01 Feb, 16:48:38
4athv38xc7s17v79ucazzr0i0/  35.39 MB  01 Feb, 16:50:56
f0bnqgf4fcr4asgi4e21tcsqa/  12.32 MB  02 Feb, 00:58:28
opqo99i4axn6fde5pjvu5hzpt/  12.32 MB  02 Feb, 00:57:35
57mr4kas99jxb9titkeackano/  11.52 MB  01 Feb, 23:07:17
mak1ji62s7axf7mespkc36g7e/  11.52 MB  01 Feb, 23:07:14
ox8718m0pwfi02zczui3qky4w/  11.52 MB  01 Feb, 23:07:16
```

And when having downloaded the data to their local machines, or simply avoiding having too much clutter on the remote storage, users can quickly delete several paths within their storage or, if they wish, remove everything. Tread carefully with the following command, but in any case, you will be asked for confirmation:

```
$ inductiva storage remove hodbisrxjxhdbkknv60xmy6ti/
Are you sure you want to remove hodbisrxjxhdbkknv60xmy6ti/? (y/[N]) y
Removing hodbisrxjxhdbkknv60xmy6ti/ in the remote storage.
Successfully removed remote path 'hodbisrxjxhdbkknv60xmy6ti/'.
```

These are the main functionalities of the CLI at the moment. In case you would like to have more access via the CLI contact us at [contact].

## What to read next
* []()
