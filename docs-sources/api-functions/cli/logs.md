# logs

The `inductiva logs` command allows you to stream the standard output (STDOUT) of a running task in real-time. This is useful for monitoring live execution progress.

## Usage

```sh
inductiva logs [-h] [--stdout] [--stderr] [--no-color] [--wait] [mode] <task_id>
```

### Description
This command will stream the logs being written to STDOUT and STDERR by a running task. 
Logs can only be streamed while the task is still in progress. If the task has finished
or transitioned to another status, logs will not be available for streaming.

## Options

- **`-h, --help`** → Show help message and exit.
- **`--stdout`** → Stream only the standard output (STDOUT) of the task.
- **`--stderr`** → Stream only the standard error (STDERR) output of the task.
- **`--no-color`** → Disable colored output for better readability in plain-text environments.
- **`--wait, -w`** → Wait for the task to start running before streaming logs. If omitted, the logs are streamed immediately or return an error if the task is not running.

## Example Usage

### Stream logs of a specific task:
```sh
inductiva logs <task_id>
```


Here's an example of how you can use the `logs` subcommand to check how your
simulation is progressing. In this case, we are trying to tracking the computation
steps, residuals, and execution times of an OpenFOAM simulation to see if everything
is evolving as expected, in realtime:

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


### Stream only the error logs of a running task:
```sh
inductiva logs <task_id> --stderr
```

### Wait for a task to start before streaming logs:
```sh
inductiva logs <task_id> --wait
```

### Stream logs without color formatting:
```sh
inductiva logs <task_id> --no-color
```

## Need Help?
Run the following command for more details:

```sh
inductiva logs --help
```

