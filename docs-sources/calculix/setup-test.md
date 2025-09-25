# Test Your Inductiva Setup
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick CalculiX simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""CalculiX example"""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "calculix-input-example.zip",
    unzip=True)

# Initialize the Simulator
calculix = inductiva.simulators.CalculiX( \
    version="2.22")

# Run simulation with config files in the input directory
task = calculix.run( \
    input_dir=input_dir,
    sim_config_filename="hueeber3.inp",
    on=cloud_machine,
    n_vcpus=4)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

<a href="https://console.inductiva.ai/editor?simulator_name=calculix" class="try-playground-button" target="_blank">
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The CalculiX simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 12/08, 09:09:28      1.081 s
	In Queue                  at 12/08, 09:09:29      58.415 s
	Preparing to Compute      at 12/08, 09:10:27      1.416 s
	In Progress               at 12/08, 09:10:29      3.056 s
		└> 2.864 s         ccx -i hueeber3
	Finalizing                at 12/08, 09:10:32      0.547 s
	Success                   at 12/08, 09:10:32      

Data:
	Size of zipped output:    216.18 KB
	Size of unzipped output:  1.87 MB
	Number of output files:   8

Estimated computation cost (US$): 0.000046 US$
```

If the task status shows **Success**, congratulations! You've successfully run a CalculiX simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

```{banner_small}
:origin: calculix_setup_test
```

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
