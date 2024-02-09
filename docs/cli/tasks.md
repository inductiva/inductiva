
### Tasks

The power of Inductiva API lies on the ability to launch multiple simulations remotely.
As an easy to track them, the CLI provides a `list` subcommand that provides information
about the last tasks launched by the user. If you wish, you can also just track
a specific task via its `--id`. Examples of usage:

```bash
$ inductiva tasks list -n 4

       ID                              SIMULATOR          STATUS         SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
       jxwt0rm8s8xspdfcegtgkkana       splishsplash       started        08 Feb, 13:25:49       08 Feb, 13:26:04       *0:00:05                 c2-standard-4
       n0zcac8rmw7xydbis3m407kb4       splishsplash       started        08 Feb, 13:25:48       08 Feb, 13:26:03       *0:00:07                 c2-standard-4
       8nmpn4h99nyfpo4da9jw2405q       splishsplash       started        08 Feb, 13:25:47       08 Feb, 13:26:02       *0:00:09                 c2-standard-4
       so6i93pi74b89rndircubp3v2       splishsplash       started        08 Feb, 13:25:47       08 Feb, 13:26:02       *0:00:10                 c2-standard-4
```

One cool tip is to use a [watch method](https://www.geeksforgeeks.org/watch-command-in-linux-with-examples/) to recurrently monitor the task information.
For example in Linux or Mac, one can monitor the task information every 10 seconds with:
```bash
$ watch -n 10 inductiva tasks list -id jxwt0rm8s8xspdfcegtgkkana
Every 10.0s: inductiva tasks list -id jxwt0rm8s8xspdfcegtgkkana                                                                                 


       ID                              SIMULATOR          STATUS         SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
       jxwt0rm8s8xspdfcegtgkkana       splishsplash       success        08 Feb, 13:25:49       08 Feb, 13:26:04       0:00:35                  c2-standard-4
```

In case, you have just a simulation and are not satisfied or you have found a mistake,
you can always `kill` it and free up the computational resources for other tasks. E.g.:

```bash
$ inductiva tasks kill cmvsc9qhz5iy86f6pef8uyxqt 
You are about to kill the following tasks:
  - cmvsc9qhz5iy86f6pef8uyxqt 
Are you sure you want to proceed (y/[N])? y
Successfully sent kill request for task cmvsc9qhz5iy86f6pef8uyxqt.
> Consider tracking the task status by calling task.get_status()
> or via the CLI using 'inductiva tasks -id cmvsc9qhz5iy86f6pef8uyxqt'.
```

Following, the instruction we confirm with:
```bash
$ inductiva tasks list -id cmvsc9qhz5iy86f6pef8uyxqt

      ID                              SIMULATOR          STATUS         SUBMITTED              STARTED         COMPUTATION TIME         RESOURCE TYPE
       cmvsc9qhz5iy86f6pef8uyxqt       splishsplash       killed         08 Feb, 13:41:06       n/a             n/a                      n/a
```

If we wanted to wait for confirmation when we sent the request to kill the task,
we can add the flag `--wait-timeout 10` in the `kill` subcommand with the number of
seconds we want to wait for confirmation.
