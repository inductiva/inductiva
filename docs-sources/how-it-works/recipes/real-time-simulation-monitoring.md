# Real-Time Simulation Monitoring with Condition-Based Termination

## The Challenge

You've just run hundreds of simulations, each generating multiple output files. Now, you need to monitor their convergence and possibly stop simulations automatically when certain criteria are met. Doing this manually is tedious, error-prone, and inefficient.

Why not automate real-time monitoring and condition-based termination instead?

## The Solution

Inductiva allows you to monitor log files in real time and trigger actions like stopping a simulation task based on output values, for example, residuals during a CFD run.

This example uses an OpenFOAM case like the one in the [OpenFOAM Quick Start guide](https://inductiva.ai/guides/openfoam/quick-start).

---

## How It Works

- The script tails the OpenFOAM log file (`log.simpleFoam`) in real time.
- It searches each line for the "Final residual" value.
- When the residual falls below a threshold (default `1e-5`), it stops (kills) the running task.

---

## Usage

Run the script passing the ID of your OpenFOAM task:

```python monitor_residuals.py <task_id>```

The script will wait for the task to start, then stream the residuals in the log file. Once the residual falls below the threshold, the task is automatically terminated.

---

## Summary

This approach helps you:

- Save computation time by terminating converged simulations automatically.
- Monitor long-running simulations in real time without manual log checking.
- Integrate with any simulation platform supported by Inductiva that writes residuals or similar convergence info to log files.

