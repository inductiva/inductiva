# Track your Tasks

The Inductiva API allows you to run multiple simulations remotely, and for every
simulation you submit to the API, you generate a `task`. The Inductiva CLI makes it
easy to monitor and manage all tasks you generated with the `list` subcommand that displays details about the most recent 
tasks you launched. You can also track a specific task using its unique`--id`.

```{seealso}
Learn more about how [tasks](https://tutorials.inductiva.ai/intro_to_api/tasks.html) are generated through the Inductiva API.
```

Here's an example of how you can use these features:

```bash
# List the last four tasks with their details
$ inductiva tasks list -n 4
       ID                              SIMULATOR          STATUS         SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
       jxwt0rm8s8xspdfcegtgkkana       splishsplash       started        08 Feb, 13:25:49       08 Feb, 13:26:04       *0:00:05                 c2-standard-4
       n0zcac8rmw7xydbis3m407kb4       splishsplash       started        08 Feb, 13:25:48       08 Feb, 13:26:03       *0:00:07                 c2-standard-4
       8nmpn4h99nyfpo4da9jw2405q       splishsplash       started        08 Feb, 13:25:47       08 Feb, 13:26:02       *0:00:09                 c2-standard-4
       so6i93pi74b89rndircubp3v2       splishsplash       started        08 Feb, 13:25:47       08 Feb, 13:26:02       *0:00:10                 c2-standard-4
```

To continuously monitor task progress, you can employ a [watch method](https://www.geeksforgeeks.org/watch-command-in-linux-with-examples/) on Linux or Mac, refreshing task information 
at set intervals.

In this example, we refreshed task information every 10 seconds:

```bash
# Monitor task status updates every 10 seconds
$ watch -n 10 inductiva tasks list --id jxwt0rm8s8xspdfcegtgkkana
Every 10.0s: inductiva tasks list --id jxwt0rm8s8xspdfcegtgkkana                                                                                 


       ID                              SIMULATOR          STATUS         SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
       jxwt0rm8s8xspdfcegtgkkana       splishsplash       success        08 Feb, 13:25:49       08 Feb, 13:26:04       0:00:35                  c2-standard-4
```
If you identify an error in your simulation or simply wish to stop it, the 
`kill` subcommand allows you to terminate it, thus freeing up the computational 
resources for other tasks:

```bash
# Terminate a specific task
$ inductiva tasks kill cmvsc9qhz5iy86f6pef8uyxqt
You are about to kill the following tasks:
  - cmvsc9qhz5iy86f6pef8uyxqt 
Are you sure you want to proceed (y/[N])? y
Successfully sent kill request for task cmvsc9qhz5iy86f6pef8uyxqt.
```
Following the termination, you can track the task status by calling task.get_status(),
or via the Inductiva CLI:

```bash
# Check the status of the terminated task
$ inductiva tasks list --id cmvsc9qhz5iy86f6pef8uyxqt

      ID                              SIMULATOR          STATUS         SUBMITTED              STARTED         COMPUTATION TIME         RESOURCE TYPE
       cmvsc9qhz5iy86f6pef8uyxqt       splishsplash       killed         08 Feb, 13:41:06       n/a             n/a                      n/a
```

If you want to wait for confirmation when you send the request to `kill` the task,
you can use the `--wait-timeout` flag with the `kill` subcommand specifying the number 
of seconds to wait for confirmation:

```bash
# Request to kill a task with a timeout for confirmation
$ inductiva tasks kill cmvsc9qhz5iy86f6pef8uyxqt --wait-timeout 10
```

Finally, the Inductiva CLI allows you to download tasks using:

```bash
$ inductiva tasks download task_1_id task_2_id
```

Additionally, we can specify specific files to download using
`--filenames` and we can specify the output directory with
`--output_dir`. When `--output_dir` is provided the tasks will be
download under `output_dir/tasks_id`.
