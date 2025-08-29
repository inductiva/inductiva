# Test Your Inductiva Setup
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick WAVEWATCH III simulation — it only takes a minute!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""WAVEWATCH III example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "wavewatch3-input-example.zip",
    unzip=True)

# Initialize the Simulator
wavewatch3 = inductiva.simulators.WaveWatch3(
     version="11-2024")

# List of commands to run
commands = ["ww3_grid", "ww3_prep", "ww3_shel"]

# Run simulation
task = wavewatch3.run(
    input_dir=input_dir,
    custom_switch="switch_PR3_UQ_MPI",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

<a href="https://console-dev.inductiva.ai/editor?simulator_name=wavewatch3" class="try-playground-button" target="_blank">
  <svg class="icon" xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
    <path d="M8 5v14l11-7z"/>
  </svg>
  Try it on our online Python Editor
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. 
The WAVEWATCH III simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 05/08, 09:07:28      0.789 s
	In Queue                  at 05/08, 09:07:29      36.552 s
	Preparing to Compute      at 05/08, 09:08:06      4.851 s
	In Progress               at 05/08, 09:08:11      70.786 s
		├> 64.133 s        bash /home/compile_ww3.sh switch_PR3_UQ_MPI true
		├> 1.071 s         ww3_grid
		├> 1.075 s         ww3_prep
		└> 4.087 s         ww3_shel
	Finalizing                at 05/08, 09:09:21      13.692 s
	Success                   at 05/08, 09:09:35      

Data:
	Size of zipped output:    636.81 MB
	Size of unzipped output:  1.48 GB
	Number of output files:   5274

Estimated computation cost (US$): 0.00060 US$
```

If the task status shows **Success**, congratulations! You've successfully run an WAVEWATCH III simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

```{banner_small}
:origin: wavewatch3
```

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
