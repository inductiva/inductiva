# Run Your First Simulation
This tutorial will show you how to run CM1 simulations using the Inductiva API. 

We will cover an example simulation available in the
[official CM1 v18 release](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r18.tar.gz),
to help you get started with simulations.

## Prerequisites

1. Download the required files
[here](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r18.tar.gz) and extract
them to a folder called `cm1r18`.

2. Create a folder named `input_files` inside the `cm1r18` directory.

Let's now copy all files we need to the `input_files` folder that you
created inside the `cm1r18` directory. Additional details about input files
can be found in the "Step 2" of this
[CM1 official tutorial](https://www2.mmm.ucar.edu/people/bryan/cm1/user_guide_brief.html).

Copy the following files to the `input_files` directory:

|            Origin           |             Description            |
|:---------------------------:|:----------------------------------:|
| `cm1r18/run/namelist.input` | Our simulation configuration file. |
| `cm1r18/run/LANDUSE.TBL` | Specifies the surface conditions, if you are using surface fluxes of heat/momentum/moisture, or if you are using the atmospheric radiation scheme. |
| `cm1r18/src/base.F` | Modifies the base-state conditions, as appropriate. There are two sections: one for the hydrostatic pressure, temperature, and moisture sounding; and one for the initial winds (u and v components). |
| `cm1r18/src/init3d.F` | Adds perturbations to the base state. Several default options are available. |
| `cm1r18/src/init_terrain.F` | Specifies the terrain via the `zs` array. |
| `cm1r18/src/init_surface.F` | Specifies the horizontal distribution of several variables, if you are using surface fluxes of heat/moisture/momentum. |

From the previous list, the only file that is **mandatory** is `namelist.input`.
The rest of the files are optional and depend on the specific configuration you
want to run. If any of the `.F` files are missing, the default files provided by
CM1 will be used.

You can also provide your own `input_sounding` file but we wont be using it in
this tutorial.

Once you have all your files, your `input_files` directory should look like this:

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

Go to https://console.inductiva.ai/tasks/4kemtaacrjjoyr92oksh819my for more details.
```

As you can see in the "In Progress" line, the part of the timeline that
represents the actual execution of the simulation, the core computation time of
this simulation was approximately 364.6 seconds (around 6 minutes).

It's that simple!

## Scaling Up Your Simulation

In order to scale your simulation to a larger machine, you need to change a
couple of lines in your `namelist.input` file and in your Python script.

Here are a list of changes you need to do:
- Change `nodex` to 4 in the `namelist.input` file.
- Change `nodey` to 4 in the `namelist.input` file.
- Change your `machine_type` to `c3d-highcpu-16` in your Python script.
- Change your `n_vcpus` to `16` in your Python script.

This is all you need to do to scale your simulation to a 16 vCPU machine.

Here are the results of running the same simulation on a few machines:

| Machine Type            | Virtual CPUs | Time             | Estimated Cost |
|-------------------------|--------------|------------------|----------------|
| **Local Ryzen 7 7700X** | 16 | 1 minute and 20 seconds | N/A |
| **Cloud c3d-highcpu-16** | 16 | 1 minute and 53 seconds | 0.0051 US$ |
| **Cloud c3d-highcpu-60** | 60 | 1 minute and 25 seconds | 0.014 US$ |

By leveraging the Inductiva API, you can efficiently scale your CM1 simulations
to meet your computational needs. Try different machine configurations and
optimize your workflow for faster, more cost-effective results!
