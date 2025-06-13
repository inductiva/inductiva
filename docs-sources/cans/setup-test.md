# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick CaNS simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""CaNS example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
	spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cansv2.4.0-input-example.zip",
    unzip=True)

# Initialize the Simulator
cans = inductiva.simulators.CaNS(\
    version="2.4.0")

# Run simulation
task = cans.run(input_dir=input_dir,
    sim_config_filename="input.nml",
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The CaNS simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 13/06, 13:35:57      0.649 s
	In Queue                  at 13/06, 13:35:58      38.191 s
	Preparing to Compute      at 13/06, 13:36:36      1.691 s
	In Progress               at 13/06, 13:36:38      4.506 s
		├> 1.184 s         mkdir -p data
		└> 3.125 s         /opt/openmpi/4.1.6/bin/mpirun --np 0 --use-hwthread-cpus cans input.nml
	Finalizing                at 13/06, 13:36:42      0.458 s
	Success                   at 13/06, 13:36:43      

Data:
	Size of zipped output:    661.22 KB
	Size of unzipped output:  2.40 MB
	Number of output files:   269

Estimated computation cost (US$): 0.000057 US$
```

If the task status shows **Success**, congratulations! You've successfully run an AMR-Wind simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

Start running simulations seamlessly! 

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
