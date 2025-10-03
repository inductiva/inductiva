# Run Your First Simulation
This tutorial will show you how to run HEC-RAS simulations using the Inductiva API. 

We will cover the `Muncie Test` from the [HEC-RAS Linux Release Notes](https://www.hec.usace.army.mil/software/hec-ras/documentation/HEC-RAS_66_Linux_Build_Release_Notes.pdf) to help you get started with simulations.

## Prerequisites
Download the required files [here](https://www.hec.usace.army.mil/software/hec-ras/downloads/Linux_RAS_v66.zip). The simulation files will be placed inside the `Linux_RAS_v66/Muncie/wrk_source` folder.

Before we begin, there's one more step: copy the
`Linux_RAS_v66/remove_HDF5_Results_Sed.py` script into the
`Linux_RAS_v66/Muncie/wrk_source` folder.

## Simulation Workflow in HEC-RAS
This simulation is executed in **six main steps**, each corresponding to a command. Below is an explanation of what each step does and why it’s necessary.

1. **Preprocess the Geometry**

   ```bash
   RasGeomPreprocess Muncie.p04.tmp.hdf x04
   ```

   This command updates the **hydraulic property tables** in the temporary geometry file (`Muncie.p04.tmp.hdf`). These tables are essential for running the unsteady flow simulation properly.

2. **Rename the Geometry File**

   ```bash
   mv Muncie.p04.tmp.hdf Muncie.p04.hdf
   ```

   The Python script in the next step always generates a file named `*.tmp.hdf`. To avoid conflicts with existing files, we rename the geometry file to `Muncie.p04.hdf`.

3. **Remove Results from the HDF5 File**

   ```bash
   python3 remove_HDF5_Results_Sed.py Muncie.p04.hdf
   ```

   This script deletes the **results data groups** (for both unsteady flow and unsteady sediment) from the HDF5 file.

   This cleanup step is necessary because there is currently **no Linux version** of the `RasUnsteadySediment` program. If these groups are not removed, the unsteady simulation will fail to start correctly.
   
4. **Run the Unsteady Flow Simulation**

   ```bash
   RasUnsteady Muncie.p04.tmp.hdf x04
   ```

   This command runs the **RAS Unsteady solver**, which performs the unsteady flow calculations based on the geometry and boundary conditions. As before, the solver outputs the results to a temporary file named `Muncie.p04.tmp.hdf`.

5. **Rename the Unsteady Output File**

   ```bash
   mv Muncie.p04.tmp.hdf Muncie.p04.hdf
   ```

   After the unsteady simulation completes, we rename the output file back to `Muncie.p04.hdf`.

6. **Run the Steady Flow Simulation**

   ```bash
   RasSteady Muncie.r04
   ```

   Finally, we run the **RAS Steady solver** using the steady simulation plan file (`Muncie.r04`). This step generates steady flow results for the same system.

## Running a HEC-RAS Simulation
Here is the code required to run a HEC-RAS simulation using the Inductiva API:

```python
"""HEC-RAS example."""
import inductiva

# Allocate Google cloud machine
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-4")

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
    input_dir="/Path/to/wrk_source",,
    commands=hec_ras_commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()
```

In this basic example, we're using a cloud machine (`c2d-highcpu-4`) equipped with 4 virtual CPUs. 
For larger or more compute-intensive simulations, consider adjusting the `machine_type` parameter to select 
a machine with more virtual CPUs and increased memory capacity. You can explore the full range of available machines [here](https://console.inductiva.ai/machine-groups/instance-types).

> **Note**: Setting `spot=True` enables the use of [spot machines](../how-it-works/machines/spot-machines.md), which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

To adapt this script for other HEC-RAS simulations, replace `input_dir` with the
path to your HEC-RAS input files.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

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

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of 
the simulation, the core computation time of this simulation was approximately 40 seconds.

```{banner_small}
:origin: hec-ras_quick_start
```