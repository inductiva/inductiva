# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. 
To confirm everything is working as expected, simply run a quick Octopus simulation — it only takes a minute!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""Octopus example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "octopus-input-example.zip",
    unzip=True)

octopus = inductiva.simulators.Octopus( \
    version="16.1")

task = octopus.run( \
    input_dir=input_dir,
    commands=["octopus"],
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The Octopus simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 14/07, 14:59:45      0.809 s
	In Queue                  at 14/07, 14:59:45      37.042 s
	Preparing to Compute      at 14/07, 15:00:23      4.96 s
	In Progress               at 14/07, 15:00:27      20.262 s
		└> 20.096 s        octopus
	Finalizing                at 14/07, 15:00:48      0.55 s
	Success                   at 14/07, 15:00:48      

Data:
	Size of zipped output:    13.57 MB
	Size of unzipped output:  21.64 MB
	Number of output files:   32

Estimated computation cost (US$): 0.00019 US$

Go to https://console.inductiva.ai/tasks/bek9bruvs0xiu1y3sa7sit3da for more details.
```

If the task status shows **Success**, congratulations! You've successfully run an Octopus simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful 
options to supercharge your simulations.

```{banner_small}
:origin: octopus
```

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!