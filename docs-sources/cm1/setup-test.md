# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick CM1 simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""CM1 example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cm1-input-example.zip",
    unzip=True)

# Initialize the Simulator
cm1 = inductiva.simulators.CM1( \
    version="21.1")

# Run simulation with config files in the input directory
task = cm1.run( \
    input_dir=input_dir,
    sim_config_filename="namelist.input",
    on=cloud_machine)

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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The CM1 simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 09/07, 11:35:56      0.813 s
	In Queue                  at 09/07, 11:35:57      36.049 s
	Preparing to Compute      at 09/07, 11:36:33      3.652 s
	In Progress               at 09/07, 11:36:36      10.312 s
		└> 10.106 s        /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus cm1.exe
	Finalizing                at 09/07, 11:36:47      0.526 s
	Success                   at 09/07, 11:36:47      

Data:
	Size of zipped output:    51.18 KB
	Size of unzipped output:  14.35 MB
	Number of output files:   16

Estimated computation cost (US$): 0.00011 US$
```

If the task status shows **Success**, congratulations! You've successfully run an CM1 simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

Start running simulations seamlessly! 

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
