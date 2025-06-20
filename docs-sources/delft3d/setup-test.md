# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick Delft3D simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""Delft3D example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "delft3d-input-example.zip",
    unzip=True)

# Initialize the Simulator
delft3d = inductiva.simulators.Delft3D(version="6.04.00")

# List of commands to run
commands = ["mpirun -np 4 d_hydro.exe config_d_hydro.xml"]

# Run simulation
task = delft3d.run( \
    input_dir=input_dir,
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

2. Open your command line, then navigate to the Desktop by running:

```
cd ~/Desktop
```

3. Execute the Python script by running:

```
python example.py
```

> **Note**: On some systems, you might need to use `python3` instead of `python`.

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The Delft3D simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal.

```
Task status: Success

Timeline:
	Waiting for Input         at 17/06, 12:33:07      0.847 s
	In Queue                  at 17/06, 12:33:08      37.008 s
	Preparing to Compute      at 17/06, 12:33:45      2.828 s
	In Progress               at 17/06, 12:33:48      6.089 s
		└> 5.941 s         mpirun -np 4 d_hydro.exe config_d_hydro.xml
	Finalizing                at 17/06, 12:33:54      0.455 s
	Success                   at 17/06, 12:33:55      

Data:
	Size of zipped output:    424.79 KB
	Size of unzipped output:  1.04 MB
	Number of output files:   15

Estimated computation cost (US$): 0.000079 US$
```

If the task status shows **Success**, congratulations! You've successfully run an Delft3D simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
