# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick AMR-Wind simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""AMR-Wind example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "amr-wind-input-example.zip",
    unzip=True)

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind( \
    version="3.4.1")

# Run simulation
task = amr_wind.run(input_dir=input_dir,
    sim_config_filename="abl_amd_wenoz.inp",
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The AMR-Wind simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 13/06, 13:26:05      0.776 s
	In Queue                  at 13/06, 13:26:06      39.502 s
	Preparing to Compute      at 13/06, 13:26:46      2.361 s
	In Progress               at 13/06, 13:26:48      4.213 s
		└> 4.074 s         /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus amr_wind abl_amd_wenoz.inp
	Finalizing                at 13/06, 13:26:52      0.663 s
	Success                   at 13/06, 13:26:53      

Data:
	Size of zipped output:    13.54 MB
	Size of unzipped output:  52.27 MB
	Number of output files:   91

Estimated computation cost (US$): 0.000065 US$
```

If the task status shows **Success**, congratulations! You've successfully run an AMR-Wind simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

Start running simulations seamlessly! 

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
