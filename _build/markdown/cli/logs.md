# logs

## inductiva logs - CLI interface

```default
inductiva logs [-h] [--stdout] [--stderr] [--no-color] [--wait] [mode]
```

The `inductiva logs` command streams the standard output (STDOUT) of a  running task in real-time, useful for monitoring live execution progress.

This command will stream the logs being written to STDOUT by a task only if the task is still in progress. If the task has finished, or transitioned to a status, the logs will not be available for streaming.

Real-time streaming of a running task’s standard error (STDERR) is also  supported via an argument.

### Positional Arguments

* [**`mode`**]() - Mode of log retrieval. Use `'SUBMITTED'` for the last submitted task, `'SUBMITTED-1'` for the  second last submitted task, and so on. `'STARTED'` and `'STARTED-n'` follow the same pattern for started tasks. Or, use a specific task ID to retrieve logs for a particular task. (default: `SUBMITTED`)

### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`--stdout`**]() - Stream the task’s standard output (STDOUT) in real time.
* [**`--stderr`**]() - Stream the task’s standard error (STDERR) in real time.
* [**`--no-color`**]() - Disable colored output.
* [**`--wait`**](), [**`-w`**]() - Wait for the task to start running before streaming logs. Without this flag, logs are streamed immediately, or an error is returned if the task is not running.

### Examples

```bash


# Real-time streaming of standard output from an OpenFOAM simulation
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
...

# Real-time streaming of a task's standard error output
$ inductiva logs f0bnqgf4fcr4asgi4e21tcsqa --stderr

# Wait for the task to start before streaming logs
$ inductiva logs f0bnqgf4fcr4asgi4e21tcsqa --wait

# Stream logs without colors
$ inductiva logs f0bnqgf4fcr4asgi4e21tcsqa --no-color
```

---
::docsbannersmall
::
