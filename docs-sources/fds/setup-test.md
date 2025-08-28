# Test Your Inductiva Setup
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick FDS simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""FDS example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-2",
    spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "fds-input-example.zip",
    unzip=True)

# Initialize the Simulator
fds = inductiva.simulators.FDS( \
    version="6.10.1")

# Run simulation
task = fds.run( \
    input_dir=input_dir,
    sim_config_filename="mccaffrey.fds",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

<a href="https://console-dev.inductiva.ai/playground?simulator_name=fds" class="try-playground-button" target="_blank">
  <svg class="icon" xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
    <path d="M8 5v14l11-7z"/>
  </svg>
  Try it in the Playground
</a>

2. Open your command line, then navigate to the Desktop by running:

```
cd ~/Desktop
```

3. Execute the Python script by running:

```
python example.py
```

> **Note**: On some systems, you might need to use `python3` instead of `python`.

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The FDS simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal.

```
Task status: Success

Timeline:
	Waiting for Input         at 21/07, 11:15:03      0.746 s
	In Queue                  at 21/07, 11:15:03      40.253 s
	Preparing to Compute      at 21/07, 11:15:44      5.44 s
	In Progress               at 21/07, 11:15:49      36.253 s
		└> 36.105 s        /opt/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux mccaffrey.fds
	Finalizing                at 21/07, 11:16:25      0.528 s
	Success                   at 21/07, 11:16:26

Data:
	Size of zipped output:    6.34 MB
	Size of unzipped output:  10.59 MB
	Number of output files:   20

Estimated computation cost (US$): 0.00017 US$
```

If the **Task status** is marked as **Success**, congratulations! You've successfully ran an FDS simulation.

You can view more details and track the full simulation progress in the [Inductiva Console](https://console.inductiva.ai/tasks).

<p align="center"><img src="./_static/set-up/console_timeline.png" alt="Task summary displayed in the Inductiva Console" width="700"></p>

This simple example tested your installation on a small machine with just 2 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

```{banner_small}
:origin: fds
```

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
