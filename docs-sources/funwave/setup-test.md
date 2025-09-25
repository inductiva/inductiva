# Test Your Inductiva Setup
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick FUNWAVE simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""FUNWAVE example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Download input files and store them in a directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "funwave-input-example.zip",
    unzip=True)

# Initialize the Simulator
funwave = inductiva.simulators.FUNWAVE( \
    version="3.6")

# Run simulation
task = funwave.run( \
    input_dir=input_dir,
    sim_config_filename="input.txt",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

<a href="https://console.inductiva.ai/editor?simulator_name=funwave" class="try-playground-button" target="_blank">
  <svg class="icon" xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
    <path d="M8 5v14l11-7z"/>
  </svg>
  Try it on our Python Editor, on any device
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The FUNWAVE simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 25/09, 10:46:34      0.887 s
	In Queue                  at 25/09, 10:46:35      33.676 s
	Preparing to Compute      at 25/09, 10:47:09      4.119 s
	In Progress               at 25/09, 10:47:13      50.9 s
		├> 1.081 s         cp /FUNWAVE-TVD-Version_3.6/Makefile .
		├> 11.084 s        make
		├> 36.113 s        /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus funwave-work/compiled_funwave input.txt
		├> 1.091 s         rm -r funwave-work
		└> 1.088 s         rm Makefile
	Finalizing                at 25/09, 10:48:04      0.555 s
	Success                   at 25/09, 10:48:04      

Data:
	Size of zipped output:    4.45 MB
	Size of unzipped output:  25.79 MB
	Number of output files:   220

Estimated computation cost (US$): 0.00036 US$
```

If the task status shows **Success**, congratulations! You've successfully run a FUNWAVE simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

```{banner_small}
:origin: funwave-setup-test
```

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!







