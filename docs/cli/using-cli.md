# CLI Usage Guide

The Inductiva CLI streamlines the management of computational resources, enabling 
you to efficiently prepare and launch the necessary environments for your simulations. 
This guide will go over how to use the Inductiva CLI for:

- [Managing your Computational Resources]()
- [Tracking your Tasks]()
- [Streaming Logs]()
- [Accessing your Remote Storage]()

## Managing your Computational Resources

The Inductiva API provides a simple way to [launch dedicated resources](../how_to/computational_resources.md) for running your simulations. With the CLI, you 
can effectively manage these computational resources, from selection and launch to 
termination, directly from your terminal. For further details on each command and 
additional options, refer to the CLI's built-in help system using the `--help` flag.

### Discovering Available Resources

Before launching any resources, you can use the CLI to go through the variety 
of machine types available and their associated costs using the `available` and 
`cost` subcommands:

```bash
$ inductiva resources available
Machine types provided in Google Cloud
```
You would get an output listing all available machines:

```{toggle}

```bash
c2: Intel Xeon Cascade Lake (2nd Gen) processor.
  > c2-standard-  [2, 4, 8, 16, 32, 60]                         

c3: Intel Xeon Sapphire Rapids (4th Gen) processor.
  > c3-highcpu-   [4, 8, 22, 44, 88, 176]                       
  > c3-standard-  [4, 8, 22, 44, 88, 176]                       
  > c3-highmem-   [4, 8, 22, 44, 88, 176]                       

h3: (Available Soon) Intel Xeon Sapphire Rapids (4th Gen) processor.
Simultaneous multithreading disabled, i.e., vCPU represents an entire core.
  > h3-standard-  [88]                                          

c2d: AMD EPYC Milan (3rd Gen) processor.
  > c2d-highcpu-  [2, 4, 8, 16, 32, 56, 112]                    
  > c2d-standard- [2, 4, 8, 16, 32, 56, 112]                    
  > c2d-highmem-  [2, 4, 8, 16, 32, 56, 112]                    

c3d: AMD EPYC Genoa (4th Gen) processor.
  > c3d-highcpu-  [4, 8, 16, 30, 60, 90, 180, 360]              
  > c3d-standard- [4, 8, 16, 30, 60, 90, 180, 360]              
  > c3d-highmem-  [4, 8, 16, 30, 60, 90, 180, 360]              

e2: Intel Xeon (up to Skylake, 1st Gen) and AMD EPYC (up to Milan, 3rd Gen) processors.
Automatically selected based on availability.
  > e2-highcpu-   [2, 4, 8, 16, 32]                             
  > e2-standard-  [2, 4, 8, 16, 32]                             
  > e2-highmem-   [2, 4, 8, 16]                                 

n2: Intel Xeon Ice Lake and Cascade Lake processors (3rd and 2nd Gen).
Cascade Lake default up to 80 vCPUs and Ice Lake for larger machines.
  > n2-highcpu-   [2, 4, 8, 16, 32, 48, 64, 80, 96]             
  > n2-standard-  [2, 4, 8, 16, 32, 48, 64, 80, 96, 128]        
  > n2-highmem-   [2, 4, 8, 16, 32, 48, 64, 80, 96, 128]        

n2d: AMD EPYC Milan or ROME processors (3rd and 2nd Gen).
  > n2d-highcpu-  [2, 4, 8, 16, 32, 48, 64, 80, 96, 128, 224]   
  > n2d-standard- [2, 4, 8, 16, 32, 48, 64, 80, 96, 128, 224]   
  > n2d-highmem-  [2, 4, 8, 16, 32, 48, 64, 80, 96]             

n1: Intel Xeon (up to Skylake, 1st Gen) processor.
Automatically selected based on availability.
  > n1-highcpu-   [1, 2, 4, 8, 16, 32, 64, 96]                  
  > n1-standard-  [1, 2, 4, 8, 16, 32, 64, 96]                  
  > n1-highmem-   [1, 2, 4, 8, 16, 32, 64, 96]  
```

If you want more detailed information you can use the `-v` flag, and you can focus
on a specific series by using the `-s` flag. 

For example:

```bash
$ inductiva resources available -s c3d -v
Machine types provided in Google Cloud
```
You'll get the following output:

```bash
c3d: AMD EPYC Genoa (4th Gen) processor.
  > c3d-highcpu-  [4, 8, 16, 30, 60, 90, 180, 360]              -> 2 GB of memory per vCPU.
  > c3d-standard- [4, 8, 16, 30, 60, 90, 180, 360]              -> 4 GB of memory per vCPU and possible local ssd integration.
  > c3d-highmem-  [4, 8, 16, 30, 60, 90, 180, 360]              -> 4 GB of memory per vCPU.
```
### Estimating Costs

You can estimate the costs of the computational resources you plan to use per hour. 
The CLI provides a cost estimation tool that considers the machine type, usage duration, 
and number of machines.

Consider the following example, where you wish to estimate the cost of four **c2-standard-8** machines:

```bash
$ inductiva resources cost c2-standard-8 --spot -n 4
```
You'll get the following output:
```bash
Estimated total cost (per machine): 0.445 (0.111) $/h.
```

### Listing Active Resources

Once you've decided and launched your resources, you can use the `list` subcommand 
to get an overview of your active computational resources:

```bash
$ inductiva resources list
```
For example, if you were running several computational resources, you'd get such
output:

```bash
Active Resources:

       NAME                                MACHINE TYPE         ELASTIC         TYPE           # MACHINES         DATA SIZE IN GB         SPOT         STARTED AT (UTC)
       api-p3kun5wyta1hacstu4xk38ujr       c2-standard-8        False           mpi            2                  10                      False        08 Feb, 12:59:10
       api-rdqprn82417bsd7id1qnac4c6       c2-standard-4        False           standard       16                 10                      False        08 Feb, 12:58:28
```

### Terminating Resources

Finally, you can terminate computational resources that are no longer needed through 
the CLI with the `terminate` subcommand. You can either choose a specific resource 
by providing its name or terminate all the resources with the `--all` flag. 
**Any of the steps require user confirmation before proceeding.** 

For example, you can choose to terminate all the resources:

```bash
$ inductiva resources terminate --all
```
In such case, you'd get the following output to confirm
your choice and terminate the resources:

```bash
You are about to terminate ALL resources.
Are you sure you want to proceed (y/[N])? y
Terminating MPICluster(name="api-p3kun5wyta1hacstu4xk38ujr"). This may take a few minutes.
MPI Cluster api-p3kun5wyta1hacstu4xk38ujr with c2-standard-8 x2 machines successfully terminated in 0:01:10.
Terminating MachineGroup(name="api-rdqprn82417bsd7id1qnac4c6"). This may take a few minutes.
Machine Group api-rdqprn82417bsd7id1qnac4c6 with c2-standard-4 machines successfully terminated in 0:01:18.
```

## Tracking your Tasks

The Inductiva API allows you to run multiple simulations remotely, and for every
simulation you submit to the API, you generate a `task`. The Inductiva CLI makes it
easy to monitor and manage all tasks you generated with the `list` subcommand that displays details about the most recent 
tasks you launched. You can also track a specific task using its unique`--id`.

:::{seealso}
Learn more about [Tasks]().
:::

Here's an example of how you can use these features:

```bash
# List the last four tasks with their details
$ inductiva tasks list -n 4
```

You'd get such output with information about the last four tasks you launched:

```bash
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
$ watch -n 10 inductiva tasks list -id jxwt0rm8s8xspdfcegtgkkana
Every 10.0s: inductiva tasks list -id jxwt0rm8s8xspdfcegtgkkana                                                                                 


       ID                              SIMULATOR          STATUS         SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
       jxwt0rm8s8xspdfcegtgkkana       splishsplash       success        08 Feb, 13:25:49       08 Feb, 13:26:04       0:00:35                  c2-standard-4
```
If you identify an error in your simulation or simply wish to stop it, the 
`kill` subcommand allows you to terminate it, thus freeing up the computational 
resources for other tasks:

```bash
# Terminate a specific task
$ inductiva tasks kill cmvsc9qhz5iy86f6pef8uyxqt
```
When you issue a `kill` command, you'll receive a prompt to confirm the action:

```bash
You are about to kill the following tasks:
  - cmvsc9qhz5iy86f6pef8uyxqt 
Are you sure you want to proceed (y/[N])? y
Successfully sent kill request for task cmvsc9qhz5iy86f6pef8uyxqt.
```
Following the termination, you can track the task status by calling task.get_status(),
or via the Inductiva CLI:

```bash
# Check the status of the terminated task
$ inductiva tasks list -id cmvsc9qhz5iy86f6pef8uyxqt
```

You'll see something like this:

```bash
      ID                              SIMULATOR          STATUS         SUBMITTED              STARTED         COMPUTATION TIME         RESOURCE TYPE
       cmvsc9qhz5iy86f6pef8uyxqt       splishsplash       killed         08 Feb, 13:41:06       n/a             n/a                      n/a
```

:::{tip}

If you want to wait for confirmation when you send the request to `kill` the task,
you can use the `--wait-timeout` flag with the `kill` subcommand specifying the number 
of seconds to wait for confirmation:

```bash
# Request to kill a task with a timeout for confirmation
$ inductiva tasks kill cmvsc9qhz5iy86f6pef8uyxqt --wait-timeout 10
```
:::
## Streaming Logs of a Task

While a task is running, one of the main concerns for scientists and engineers is
making sure the simulation is progressing as expected or identifying any potential 
errors. To help with this, the Inductiva CLI provides a `logs` subcommand that enables
you to stream the logs of a running task in **real time**.

In this way, running a simulation on a powerful remote machine will feel exactly 
like running it on your local computer!

Here's an example of how you can use the `logs` subcommand to stay updated on your simulation's 
progress:

```bash
# Stream logs in real time for a specific task
$ inductiva logs f0bnqgf4fcr4asgi4e21tcsqa
```

This command connects you directly to the ongoing logs of the specified task, displaying updates such as computation steps, residuals, and execution times as 
if you were running the simulation on your own machine:


```bash
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
Combining both subcommand `log` and `kill` will allow you to quickly adjust or 
terminate simulations to achieve the desired simulation configurations.

At the end of the simulation, the log stream will pause without new messages. In 
this case, you can exit the log stream by simply using `Ctrl+C`, which stops the 
log monitoring without affecting the simulation or its outputs.

## Accessing your Remote Storage

After your simulations finish, the results are securely stored in your **remote 
storage bucket**. The Inductiva CLI connects you to your remote bucket, and
provides you with various commands to effectively manage all the data of your 
simulation outputs.

### Checking Storage Usage

You can check how much storage space your simulations are currently occupying:

```bash
$ inductiva storage size
```
You'll get such output:

```bash
Total user's remote storage in use: 2.79 GB
```
### Listing Storage Contents

You can have a detailed view of what's in your storage, including sorting options 
for easier navigation:

```bash
$ inductiva storage list --max-results 10 --order-by size --sort-order desc
```
You'll get such overview that lists the contents of your storage, allowing you to 
see the largest files or the most recent ones, depending on your sorting preferences:

```bash
       NAME                             SIZE           CREATION TIME
       0bet8jrpp2gz974n42nsd9n2p/       56.11 MB       06 Feb, 11:32:29
       05ujj5m0ytdkckxwk1tq1b5io/       27.93 MB       08 Feb, 09:19:44
       6a2h1wnxywpea8jfoxiikdjf7/       26.49 MB       07 Feb, 13:47:03
       f8joznwc9xf9a4nypcaei6v2s/       12.79 MB       07 Feb, 09:16:55
       dpq2cv6b5f9p1c77nc8anjo10/       12.00 MB       08 Feb, 09:39:31
       r4kerxf4b53krgn0s3fyece3b/       11.92 MB       07 Feb, 11:47:48
       j9qzrpiohgt7x97od3tw4wccd/       11.74 MB       07 Feb, 11:47:46
       iqi71gonoacfj7fknox3rvnq2/       11.52 MB       07 Feb, 11:47:45
       dxmnxdrfrv84pfbzbvm9v0dat/       11.43 MB       07 Feb, 11:47:43
       bgtwgnnyq5qa5hecegzdx6okr/       11.36 MB       07 Feb, 11:47:40
```

### Cleaning Up Storage
Once you've backed up your data locally or need to free up space, you can delete 
data from your remote storage. You can delete several paths within your
storage, or even delete everything with the `remove` command.

Here's an example where you remove a path within your storage: 

```bash
$ inductiva storage remove hodbisrxjxhdbkknv60xmy6ti/
```
Use this command with caution as it permanently deletes data from your remote 
storage. **The CLI always prompts for confirmation to prevent accidental data loss**:

```bash
You are about to remove the following paths from your remote storage space:
  - 0bet8jrpp2gz974n42nsd9n2p/
Are you sure you want to proceed (y/[N])? y
Removing 0bet8jrpp2gz974n42nsd9n2p/ in the user's remote storage.
Successfully removed remote path '0bet8jrpp2gz974n42nsd9n2p/'.
```
You can alternatively clear all storage by adding the `--all` flag to the `remove` command.
Such commands will always be followed with a confirmation prompt to ensure user intention
and prevent irreversible loss.

These are the main functionalities of the CLI for storage at the moment. In case
you would like to have more functionalities via the CLI [contact us](support@inductiva.ai).
