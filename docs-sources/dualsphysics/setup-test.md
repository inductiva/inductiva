# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick DualSPHysics simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""DualSPHysics example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "dualsphysics-input-example.zip",
    unzip=True)

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics( \
    version="5.2.1")

# Run simulation 
task = dualsphysics.run( \
    input_dir=input_dir,
    shell_script="run.sh",
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The DualSPHysics simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal.

```
Task status: Success

Timeline:
	Waiting for Input         at 17/06, 10:39:13      0.678 s
	In Queue                  at 17/06, 10:39:13      38.199 s
	Preparing to Compute      at 17/06, 10:39:51      1.999 s
	In Progress               at 17/06, 10:39:53      227.585 s
		└> 227.43 s        bash run.sh
	Finalizing                at 17/06, 10:43:41      1.508 s
	Success                   at 17/06, 10:43:43      

Data:
	Size of zipped output:    64.81 MB
	Size of unzipped output:  117.44 MB
	Number of output files:   265

Estimated computation cost (US$): 0.0017 US$
```

If the task status shows **Success**, congratulations! You've successfully run an DualSPHysics simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
