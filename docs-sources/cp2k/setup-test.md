# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick CP2K simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""CP2K example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cp2k-input-example.zip",
    unzip=True)

# Initialize the Simulator
cp2k = inductiva.simulators.CP2K( \
    version="2025.1")

# Run simulation
task = cp2k.run( \
    input_dir=input_dir,
    sim_config_filename="Ac.inp",
    use_hwthread=True,
    n_vcpus=4,
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The CP2K simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal.

```
Task status: Success

Timeline:
	Waiting for Input         at 17/06, 11:34:32      0.776 s
	In Queue                  at 17/06, 11:34:32      38.612 s
	Preparing to Compute      at 17/06, 11:35:11      10.197 s
	In Progress               at 17/06, 11:35:21      5.224 s
		└> 5.072 s         /opt/openmpi/5.0.6/bin/mpirun --np 4 --use-hwthread-cpus cp2k.psmp Ac.inp
	Finalizing                at 17/06, 11:35:26      0.369 s
	Success                   at 17/06, 11:35:27      

Data:
	Size of zipped output:    7.05 KB
	Size of unzipped output:  37.39 KB
	Number of output files:   2

Estimated computation cost (US$): 0.00012 US$
```

If the task status shows **Success**, congratulations! You've successfully run an CP2K simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
