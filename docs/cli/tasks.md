
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