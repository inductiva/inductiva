# Real-Time Simulation Monitoring with Condition-Based Termination

## The Challenge

You've just run hundreds of simulations, each generating multiple output files. Now, you need to monitor their convergence and possibly stop simulations automatically when certain criteria are met. Doing this manually is tedious, error-prone, and inefficient.

Why not automate real-time monitoring and condition-based termination instead?

## The Solution

Inductiva allows you to monitor log files in real time and trigger actions like stopping a simulation task based on output values, for example, residuals during a CFD run.

Below is a Python script example for monitoring the residuals of an OpenFOAM simulation. When the residual value drops below a given threshold, the script automatically kills the task to save computing resources.

This example uses an OpenFOAM case like the one in the [OpenFOAM Quick Start guide](https://inductiva.ai/guides/openfoam/quick-start).

## Script monitor_residuals.py

```python
import re
import argparse
import inductiva


class ResidualFilter:

    def __init__(self, pattern, threshold, task):
        self.pattern = pattern
        self.threshold = threshold
        self.task = task
        self.stop = False

    def write(self, text):
        if self.stop:
            return

        for line in text.splitlines():
            match = self.pattern.search(line)
            if match:
                residual_value = float(match.group(1))
                print(f"Residual found: {residual_value}")
                if residual_value < self.threshold:
                    self.stop = True
                    print("Residual below threshold, killing task...")
                    self.task.kill()

    def flush(self):
        pass  # For compatibility


def main():
    parser = argparse.ArgumentParser(
        description="Monitor and kill a task based on residual values.")
    parser.add_argument("task_id", help="ID of the task.")
    parser.add_argument("--filename",
                        default="log.simpleFoam",
                        help="Log file to monitor (default: log.simpleFoam).")
    parser.add_argument("--threshold",
                        type=float,
                        default=1e-5,
                        help="Residual threshold to trigger task kill.")

    args = parser.parse_args()

    pattern = re.compile(r"Final residual\s*=\s*([\d.eE+-]+)")
    task = inductiva.tasks.Task(args.task_id)

    task.wait_for_status("computation-started")
    print(f"Tailing {args.filename} for residuals...")

    task.tail_files(
        tail_files=[args.filename],
        lines=10,
        follow=True,
        wait=True,
        fout=ResidualFilter(pattern, args.threshold, task),
    )

if __name__ == "__main__":
    main()
```


## How It Works

- The script tails the OpenFOAM log file (`log.simpleFoam`) in real time.
- It searches each line for the "Final residual" value.
- When the residual falls below a threshold (default `1e-5`), it stops (kills) the running task.


## Usage

Run the script passing the ID of your OpenFOAM task:

```python monitor_residuals.py <task_id>```

The script will wait for the task to start, then stream the residuals in the log file. Once the residual falls below the threshold, the task is automatically terminated.


## Summary

This approach helps you:

- Save computation time by terminating converged simulations automatically.
- Monitor long-running simulations in real time without manual log checking.
- Integrate with any simulation platform supported by Inductiva that writes residuals or similar convergence info to log files.

