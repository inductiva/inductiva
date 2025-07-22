# Real-Time Monitoring & Conditional Auto Termination

You’re running one or more simulations that take time and resources, and you want to stop them automatically once they converge based on a condition (e.g., residuals falling below a threshold). Doing this manually is tedious, error-prone, and inefficient.

Why not automate real-time monitoring and condition-based termination instead?

## The Solution

This recipe allows **real-time log monitoring** and **automatic stopping** of simulations when a condition is met, such as convergence. It applies to **any simulator** that writes convergence-related values (e.g., residuals) to a log file in a consistent format.

Benefits:
- Save computation time by terminating converged simulations automatically.
- Monitor long-running simulations in real time without manual log checking.
- Works with any simulator that logs convergence info, like residuals, to a file.


## Requirements:
- Your simulator should **write convergence indicators** (e.g., residuals) to a **log file** during execution.
- The relevant values must be parseable using a regex pattern (can be adjusted in the script).
- You must know the `task_id` of the running simulation.

This example is designed for the [OpenFOAM Quick Start guide](https://inductiva.ai/guides/openfoam/quick-start), but the approach is general.


## How to Use It in Your Workflow
First, launch your Inductiva simulation. Once you have the `task_id` (you can get it from the Python API or copy it from the web console), you’ll need to run the monitoring script (shown below) **in parallel** to start tracking the log file.

```
python my_simulation.py  # Your script that launches the simulation
python monitor_residuals.py <task_id>  # This script monitors the simulation
```

> 💡 Make sure to run the monitoring script while the simulation is still running, so it can tail the log file in real time.


This script:
- Tails the simulation log file (`log.simpleFoam`) in real time.
- Searches each line for a residual value using a regular expression.
- If the residual falls below a threshold (default: `1e-5`), it automatically stops the running task to save resources.


## Python Script monitor_residuals.py
The following script can be used to monitor residuals during a simulation run. When the residual value drops below a given threshold, it automatically stops the task to save computing resources.

Copy and save it to a file named `monitor_residuals.py`.


```python
import re
import argparse
import inductiva


class ResidualFilter:

    def __init__(self, pattern, threshold, task):
        self.pattern = pattern          # Regex pattern to extract residuals
        self.threshold = threshold      # Residual threshold to stop the task
        self.task = task                # The Inductiva task being monitored
        self.stop = False               # Flag to avoid stopping more than once

    def write(self, text):
        """Function called with new log lines."""

        if self.stop:
            return

        for line in text.splitlines():
            match = self.pattern.search(line) # Search for regex pattern
            if match:
                residual = float(match.group(1))
                print(f"Residual found: {residual}")
                if residual < self.threshold: # Check if residual below threshold
                    self.stop = True
                    print("Residual below threshold, stopping task...")
                    self.task.kill() # Stop the task

    def flush(self):
        pass  # For compatibility


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Monitor and stop a task based on residual values.")
    parser.add_argument("task_id", help="ID of the task.")
    parser.add_argument("--filename",
                        default="log.simpleFoam",
                        help="Log file to monitor (default: log.simpleFoam).")
    parser.add_argument("--threshold",
                        type=float,
                        default=1e-5,
                        help="Residual threshold to trigger task stop.")

    args = parser.parse_args()

    # Compile the regex to extract "Final residual = ..." from log lines
    pattern = re.compile(r"Final residual\s*=\s*([\d.eE+-]+)")
    
    task = inductiva.tasks.Task(args.task_id) # Load the task using its ID

    # Wait until the task actually starts computation
    task.wait_for_status("computation-started")
    print(f"Tailing {args.filename} for residuals...")

    # Stream the log file and pipe output to our ResidualFilter
    task.tail_files(
        tail_files=[args.filename],
        lines=10, # Start tailing from the last 10 lines
        follow=True, # Keep reading new lines in real time
        wait=True,  # Wait for the log file to be available
        fout=ResidualFilter(pattern, args.threshold, task),
    )

if __name__ == "__main__":
    main()
```

> 💡 If your simulator uses a different log format, adjust the regular expression in the script accordingly.
