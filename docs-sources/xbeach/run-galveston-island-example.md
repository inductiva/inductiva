# Run the Galveston Island Beach and Dune Simulation
In this tutorial we will show you how to use the Inductiva API to run an advanced XBeach case that requires significant computing power.

## Objective
The goal of this tutorial is to demonstrate how to run the `Galveston Island` use case from the [GRIIDC repository](https://www.griidc.org/), a repository at Texas A&M University-Corpus Christi’s Harte Research Institute for Gulf of Mexico Studies.

## Prerequisites
1. We will be using the dataset titled "XBeach model setup and results for beach and dune enhancement scenarios on Galveston Island, Texas," which is available [here](https://data.griidc.org/data/HI.x833.000:0001). Visit the [dataset page](https://data.griidc.org/data/HI.x833.000:0001#individual-files).
2. Navigate to **Files >> XBeach_Model_Runs >> Beach_Nourish_Only >> Input** and download all the files in this directory.
3. Save these files in a folder named `Beach_Nourish_Only`, ensuring the directory structure is as follows:

```
	ls -las Beach_Nourish_Only 
	total 130976
    0 drwxr-xr-x   12 paulobarbosa  staff       384 Nov  6 10:15 .
    0 drwx------@ 124 paulobarbosa  staff      3968 Nov  6 10:14 ..
    8 -rw-r--r--@   1 paulobarbosa  staff      3069 Nov  6 10:14 README.txt
26184 -rw-r--r--@   1 paulobarbosa  staff  13404906 Nov  6 10:14 bed.dep
26184 -rw-r--r--@   1 paulobarbosa  staff  13404906 Nov  6 10:14 bedfricfile.txt
   16 -rw-r--r--@   1 paulobarbosa  staff      5324 Nov  6 10:13 jonswap3.txt
26184 -rw-r--r--@   1 paulobarbosa  staff  13404906 Nov  6 10:14 nebed.dep
    8 -rw-r--r--@   1 paulobarbosa  staff      2296 Nov  6 10:15 params.txt
   16 -rw-r--r--@   1 paulobarbosa  staff      4850 Nov  6 10:14 tide.txt
26184 -rw-r--r--@   1 paulobarbosa  staff  13404906 Nov  6 10:14 x.grd
    8 -rw-r--r--@   1 paulobarbosa  staff       635 Nov  6 10:14 xbeach.slurm
26184 -rw-r--r--@   1 paulobarbosa  staff  13404906 Nov  6 10:14 y.grd
```

## Adjust Simulation Parameters
For a faster simulation, modify the `params.txt` file as follows:
- Add `single_dir = 0` after the header (required for XBeach v10+).
- Set `tstop` to `34560` to shorten the simulation duration.

## Run Your Simulation
Here is the code required to run an XBeach simulation using the Inductiva API:

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
	provider="GCP",
    machine_type="c3d-highcpu-90",
    spot=True,
    data_disk_gb=20)
cloud_machine.start()

# Initialize the Simulator
XBeach = inductiva.simulators.XBeach( \
    version="1.24")

# Run simulation 
task = xbeach.run(
    input_dir="Beach_Nourish_Only",
    sim_config_filename="params.txt",
    n_vcpus=90,
    on=cloud_machine)

# task.wait() is a blocking call and will only return once the simulation
# ends. However, you can close your terminal without interrupting the 
# simulation. To check the status of the simulation, you can use Inductiva # CLI (Command Line Interface) tools from another terminal.
task.wait()

# Terminate your allocated cloud machine at the end of the simulation
cloud_machine.terminate()

# Let's get a small summary of the run
task.print_summary()
```

This simulation runs on a `c3d-highcpu-90` machine with a 20 GB disk. 

**Note**: `spot` machines are a lot cheaper but may be terminated by the
	provider if necessary.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 10/04, 08:42:18      2.602 s
	In Queue                  at 10/04, 08:42:21      67.241 s
	Preparing to Compute      at 10/04, 08:43:28      1.616 s
	In Progress               at 10/04, 08:43:30      3066.09 s
		└> 3065.902 s      /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 90 xbeach params.txt
	Finalizing                at 10/04, 09:34:36      14.752 s
	Success                   at 10/04, 09:34:51      

Data:
	Size of zipped output:    405.02 MB
	Size of unzipped output:  668.23 MB
	Number of output files:   28

Estimated computation cost (US$): 0.73 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, the core computation time 
of this simulation was approximately 51 minutes.

It's that simple! 🚀




