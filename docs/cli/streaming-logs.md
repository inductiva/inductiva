# Stream Logs of a Task

While a task is running, one of the main concerns for scientists and engineers is
making sure the simulation is progressing as expected or identifying any potential 
errors. To help with this, the Inductiva CLI provides a `logs` subcommand that enables
you to stream the logs of a running task in **real time**.

In this way, running a simulation on a powerful remote machine will feel exactly 
like running it on your local computer!

## Check the Log Stream

Here's an example of how you can use the `logs` subcommand to stay updated on your simulation's 
progress. This command connects you directly to the ongoing logs of the specified 
task, displaying updates such as computation steps, residuals, and execution times as 
if you were running the simulation on your own machine:

```bash
# Stream logs in real time for a specific task
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
Combining both subcommand `log` and `kill` will allow you to quickly adjust or 
terminate simulations to achieve the desired simulation configurations.

## Save the Log Stream

When the simulation ends, the log stream will automatically close, indicating
that there are no new messages. The stream consumer's shutdown doesn't impact the 
simulation or its outputs, and you don't need to take any further action.

You can save the log stream content and the stream consumer's status for later 
inspection by redirecting them to a file:

- **Redirect stdout to a file**: You can redirect stdout (file descriptor 1) to 
a file with the redirection operator (`>`) to save the entire log stream content. 
For example:

    ```bash
    $ inductiva logs TASK_ID 1>out.txt
    ```
- **Redirect stderr to a file**: You can redirect the status of the stream consumer, 
outputted to stderr (file descriptor 2), to a file like this:

    ```bash
    $ inductiva logs TASK_ID 2>err.txt
    ```
- **Redirect stdout and stderr to separate files**: You can redirect both stdout 
and stderr simultaneously to separate files using:

    ```bash
    $ inductiva logs TASK_ID 1>out.txt 2>err.txt
    ```
- **Disable ANSI globally**: If you need to disable ANSI escape codes, which are 
used to support the status bar at the bottom of the log stream, either export the 
`ANSI_ENABLED` environment variable or set it locally:

    ```bash
    $ export ANSI_ENABLED=0 # Applies to the entire shell session
    $ ANSI_ENABLED=0 inductiva logs TASK_ID. # Applies only to this command
    ```
