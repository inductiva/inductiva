# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick XBeach simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""XBeach example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "xbeach-input-example.zip",
    unzip=True)

# Initialize the Simulator
xbeach = inductiva.simulators.XBeach( \
    version="1.24")

# Run simulation
task = xbeach.run( \
    input_dir=input_dir,
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The XBeach simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 27/06, 18:45:04      0.792 s
	In Queue                  at 27/06, 18:45:05      38.483 s
	Preparing to Compute      at 27/06, 18:45:43      3.846 s
	In Progress               at 27/06, 18:45:47      47.247 s
		└> 47.114 s        /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus xbeach params.txt
	Finalizing                at 27/06, 18:46:34      2.599 s
	Success                   at 27/06, 18:46:37      

Data:
	Size of zipped output:    74.60 MB
	Size of unzipped output:  82.73 MB
	Number of output files:   12

Estimated computation cost (US$): 0.00039 US$
```

If the task status shows **Success**, congratulations! You've successfully run an XBeach simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

Start running simulations seamlessly! 

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!