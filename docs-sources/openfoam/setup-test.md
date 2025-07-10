# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. 
To confirm everything is working as expected, simply run a quick OpenFOAM simulation — it only takes a minute!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""OpenFOAM example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "openfoam-esi-input-example.zip",
    unzip=True)

# Initialize the Simulator
openfoam = inductiva.simulators.OpenFOAM( \
    distribution="esi",
    version="2406")

# Run simulation 
task = openfoam.run( \
    input_dir=input_dir,
    shell_script="./Allrun",
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The OpenFOAM simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 10/07, 00:02:39      1.973 s
	In Queue                  at 10/07, 00:02:41      33.46 s
	Preparing to Compute      at 10/07, 00:03:14      4.623 s
	In Progress               at 10/07, 00:03:19      64.347 s
		└> 64.188 s        bash ./Allrun
	Finalizing                at 10/07, 00:04:23      1.999 s
	Success                   at 10/07, 00:04:25      

Data:
	Size of zipped output:    99.56 MB
	Size of unzipped output:  155.56 MB
	Number of output files:   223

Estimated computation cost (US$): 0.00052 US$
```

If the task status shows **Success**, congratulations! You've successfully run an OpenFOAM simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful 
options to supercharge your simulations.

Start running simulations seamlessly! 

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!