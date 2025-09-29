# Test Your Inductiva Setup
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick Elmer simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""Elmer example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "elmer-input-example.zip",
    unzip=True)

# Initialize the Simulator
elmer = inductiva.simulators.Elmer( \
    version="9.0")

# Run simulation with config files in the input directory
task = elmer.run( \
    input_dir=input_dir,
    commands=[
        "ElmerGrid 1 2 winkel.grd",
        "ElmerSolver case.sif -ipar 2 1 1"],
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()

```

<a href="https://console.inductiva.ai/editor?simulator_name=elmer" class="try-playground-button" target="_blank">
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The Elmer simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 19/09, 11:38:01      0.699 s
	In Queue                  at 19/09, 11:38:02      36.602 s
	Preparing to Compute      at 19/09, 11:38:39      4.037 s
	In Progress               at 19/09, 11:38:43      14.373 s
		├> 1.073 s         ElmerGrid 1 2 winkel.grd
		└> 13.081 s        ElmerSolver case.sif -ipar 2 1 1
	Finalizing                at 19/09, 11:38:57      0.534 s
	Success                   at 19/09, 11:38:58      

Data:
	Size of zipped output:    2.71 MB
	Size of unzipped output:  6.25 MB
	Number of output files:   11

Estimated computation cost (US$): 0.00014 US$
```

If the task status shows **Success**, congratulations! You've successfully run a Elmer simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

```{banner_small}
:origin: elmer_setup_test
```

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
