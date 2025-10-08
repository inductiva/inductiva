# Test Your Inductiva Setup
Before diving into tutorials and benchmarks, let's ensure that your Inductiva Python package is properly set up. To confirm everything is working as expected, simply run a quick HEC-RAS simulation — it only takes a few seconds!

## Step 1: Copy and Run the Code

1. Copy the code below and save it as `example.py` on your Desktop (or in your preferred directory).

```python
"""HEC-RAS example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "hec-ras-input-example.zip",
    unzip=True)

# Initialize the Simulator
hec_ras = inductiva.simulators.HecRas( \
    version="6.6")

# Specify the HEC-RAS commands you want to run, separated by commas
hec_ras_commands = [
        'RasGeomPreprocess Muncie.p04.tmp.hdf x04',
        'mv Muncie.p04.tmp.hdf Muncie.p04.hdf',
        'python3 remove_HDF5_Results_Sed.py Muncie.p04.hdf',
        'RasUnsteady Muncie.p04.tmp.hdf x04',
        'mv Muncie.p04.tmp.hdf Muncie.p04.hdf',
        'RasSteady Muncie.r04']

# Run simulation
task = hec_ras.run( \
    input_dir=input_dir,
    commands=hec_ras_commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
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

All the necessary simulation artifacts and configuration files will be automatically downloaded to your computer. The HEC-RAS simulation will then be sent to a cloud machine for execution.

## Step 2: Verify the Task Status
After the simulation completes, a task summary will be displayed in your terminal, as shown below. 

```
Task status: Success

Timeline:
	Waiting for Input         at 15/09, 11:13:11      1.244 s
	In Queue                  at 15/09, 11:13:12      33.81 s
	Preparing to Compute      at 15/09, 11:13:46      4.372 s
	In Progress               at 15/09, 11:13:50      39.885 s
		├> 1.928 s         RasGeomPreprocess Muncie.p04.tmp.hdf x04
		├> 1.082 s         mv Muncie.p04.tmp.hdf Muncie.p04.hdf
		├> 1.078 s         python3 remove_HDF5_Results_Sed.py Muncie.p04.hdf
		├> 31.105 s        RasUnsteady Muncie.p04.tmp.hdf x04
		├> 1.09 s          mv Muncie.p04.tmp.hdf Muncie.p04.hdf
		└> 3.082 s         RasSteady Muncie.r04
	Finalizing                at 15/09, 11:14:30      0.839 s
	Success                   at 15/09, 11:14:31      

Data:
	Size of zipped output:    24.52 MB
	Size of unzipped output:  36.50 MB
	Number of output files:   9

Estimated Task Compute Cost = 0.00032 US$
Task Orchestration Fee = 0.01 US$
Total Estimated Cost = 0.01032 US$
Learn more about costs at: https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost
```

If the task status shows **Success**, congratulations! You've successfully run a HEC-RAS simulation.

This simple example tested your installation on a small machine with just 4 virtual CPUs. Inductiva offers far more powerful options to supercharge your simulations.

```{banner_small}
:origin: hec-ras_setup_test
```

## Need Help?
If you encounter any issues or need further assistance, don't hesitate to [**Contact Us**](mailto:support@inductiva.ai). We're here to help!
