# Run Your First Simulation
This tutorial will show you how to run CM1 simulations using the Inductiva API. 

We will cover an example simulation available in the
[official CM1 v18 release](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r18.tar.gz),
to help you get started with simulations.

## Prerequisites
Before running the simulation, you need to prepare the required files and directory structure.

1. **Download the Files**
Download the files [here](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r18.tar.gz).

2. **Create the Input Directory**
Inside the `cm1r18` folder, create a new subfolder named `input_files`.

3. **Copy the Required Files**
Now, copy the following files into the `input_files` directory. For more information about these input files, refer to **Step 2** 
of the [CM1 official guide](https://www2.mmm.ucar.edu/people/bryan/cm1/user_guide_brief.html).

| File Path                    | Description                                                                |
|:----------------------------:|:--------------------------------------------------------------------------:|
| `cm1r18/run/namelist.input`  | Simulation configuration file. **Mandatory**.                          |
| `cm1r18/run/LANDUSE.TBL`     | Defines surface conditions used in radiation and flux schemes.             |
| `cm1r18/run/LANDUSE.TBL`     | Defines surface conditions used when enabling surface fluxes (heat, momentum, moisture) or the atmospheric radiation scheme.                        |
| `cm1r18/src/base.F`          | Configures base-state conditions. Includes two parts: (1) hydrostatic pressure, temperature, and moisture sounding, and (2) initial winds (u, v). |
| `cm1r18/src/init_terrain.F`  | Specifies terrain via the `zs` array.                                      |
| `cm1r18/src/init_surface.F`  | Sets horizontal surface distributions (e.g., heat, moisture).              |

Of the files listed above, only `namelist.input` is required.
The other files are optional and depend on the specific configuration of your simulation. If any of the `.F` files are not provided, the default files provided by CM1 will be used.

You may also include a custom `input_sounding` file, although we won’t be using one in this tutorial.

Once all necessary files are in place, your `input_files` directory should look like this:

```
-rw-r--r--@ 1 paulobarbosa  staff   5125 Jul 26  2015 LANDUSE.TBL
-rw-r--r--@ 1 paulobarbosa  staff  63930 Oct  7  2015 base.F
-rw-r--r--@ 1 paulobarbosa  staff  53493 Oct  7  2015 init3d.F
-rw-r--r--@ 1 paulobarbosa  staff   7952 Aug 31  2015 init_surface.F
-rw-r--r--@ 1 paulobarbosa  staff   8952 Aug 13  2015 init_terrain.F
-rw-r--r--@ 1 paulobarbosa  staff   6193 Oct  6  2015 namelist.input
```

Since this is just a introduction tutorial we won't be editing any of the files
and will use the default configuration.

> Note: In this case, because we are not editing any of the files we could
simply just send the file `namelist.input` in the `input_files` directory. We
added all the other files as an example of what you could do.

## Running a CM1 Simulation
Here is the code required to run a CM1 simulation using the Inductiva API:

```python
"""CM1 example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
	spot=True)

# Initialize the Simulator
cm1 = inductiva.simulators.CM1(
    version="18")

# Run simulation with config files in the input directory
task = cm1.run(input_dir="/Path/to/input_files",
               sim_config_filename="namelist.input",
			   # optional config files
               base="base.F",
               init3d="init3d.F",
               init_surface="init_surface.F",
               init_terrain="init_terrain.F",
               landuse="LANDUSE.TBL",
               n_vcpus=1,
               on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the
provider if necessary.

To adapt the code for this or any other use case, simply replace `input_dir`
with the path to your CM1 input files and set the `sim_config_filename`
accordingly.

When the simulation is complete, we terminate the machine, download the results
and print a summary of the simulation as shown below.

```
inductiva tasks info 4kemtaacrjjoyr92oksh819my

Task status: Success

Timeline:
	Waiting for Input         at 13/03, 15:00:21      0.918 s
	In Queue                  at 13/03, 15:00:22      17.618 s
	Preparing to Compute      at 13/03, 15:00:40      1.823 s
	In Progress               at 13/03, 15:00:42      364.628 s
		├> 1.164 s         cp -r /cm1 /workdir/output/artifacts/__cm1
		├> 1.151 s         cp -f base.F /workdir/output/artifacts/__cm1/src
		├> 1.151 s         cp -f init3d.F /workdir/output/artifacts/__cm1/src
		├> 1.081 s         cp -f init_terrain.F /workdir/output/artifacts/__cm1/src
		├> 1.066 s         cp -f init_surface.F /workdir/output/artifacts/__cm1/src
		├> 59.115 s        make -C /workdir/output/artifacts/__cm1/src
		├> 298.356 s       /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 1 /workdir/output/artifacts/__cm1/run/cm1.exe namelist.input
		└> 1.065 s         rm -r /workdir/output/artifacts/__cm1
	Finalizing                at 13/03, 15:06:46      2.026 s
	Success                   at 13/03, 15:06:48      

Data:
	Size of zipped output:    39.11 MB
	Size of unzipped output:  128.26 MB
	Number of output files:   13

Estimated computation cost (US$): 0.0043 US$
```

As you can see in the "In Progress" line, the part of the timeline that
represents the actual execution of the simulation, the core computation time of
this simulation was approximately 364.6 seconds (around 6 minutes).

It's that simple!

## Scaling Up Your Simulation
To run your simulation on a larger machine, you’ll need to make a few small changes to both your `namelist.input` file 
and your Python script.

### Required Changes
Update the following parameters:

* In `namelist.input`:
	- Set `nodex` = 4
	- Set `nodey` = 4
* In your Python script:
	- Set `machine_type` = "c3d-highcpu-16"
	- Set `n_vcpus` = 16

That’s all it takes to scale your simulation to a 16 vCPU machine.

### Performance Comparison
Here are the results of running the same simulation on different machines:

| Machine Type             | Virtual CPUs     | Time             | Estimated Cost (USD) |
|--------------------------|------------------|------------------|----------------------|
| **Local Ryzen 7 7700X**  | 16               | 1 min, 20s       | N/A                  |
| **Cloud c3d-highcpu-16** | 16               | 1 min, 3s        | 0.0051               |
| **Cloud c3d-highcpu-60** | 60               | 1 min, 25s       | 0.014                |

With the **Inductiva API**, you can easily scale your CM1 simulations to match your computational demands. Whether you need faster runtimes or lower costs, experimenting with different machine configurations allows you to find the optimal balance for your workflow.