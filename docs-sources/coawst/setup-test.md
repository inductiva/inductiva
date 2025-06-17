# Test Your Inductiva Setup ⚙️
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick COAWST simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""COAWST example."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Download the input files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "coawst-input-example.zip",
    unzip=True)

# Initialize the Simulator
coawst = inductiva.simulators.COAWST( \
    version="3.8")

# Run simulation
task = coawst.run( \
    input_dir=input_dir,
    sim_config_filename="ocean_ducknc.in",
    build_coawst_script="build_coawst.sh",
    n_vcpus=4,
    use_hwthread=True,
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The COAWST simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal.

```
Task status: Success

Timeline:
	Waiting for Input         at 17/06, 10:30:09      0.966 s
	In Queue                  at 17/06, 10:30:10      35.466 s
	Preparing to Compute      at 17/06, 10:30:46      10.12 s
	In Progress               at 17/06, 10:30:56      310.223 s
		├> 18.086 s        cp -r /opt/COAWST /workdir/output/artifacts/__COAWST
		├> 1.077 s         create_all_sim_links
		├> 120.181 s       bash build_coawst.sh
		├> 168.245 s       /opt/openmpi/4.1.6/bin/mpirun --np 4 --use-hwthread-cpus coawstM ocean_ducknc.in
		├> 1.072 s         rm -r __COAWST
		└> 1.075 s         clean_all_sim_links
	Finalizing                at 17/06, 10:36:06      0.89 s
	Success                   at 17/06, 10:36:07      

Data:
	Size of zipped output:    12.23 MB
	Size of unzipped output:  133.15 MB
	Number of output files:   8

Estimated computation cost (US$): 0.0023 US$
```

If the task status shows **Success**, congratulations! You've successfully run an COAWST simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
