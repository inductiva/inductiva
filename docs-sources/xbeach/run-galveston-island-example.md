# Run the Galveston Island Beach and Dune Simulation
This tutorial walks you through running a high-fidelity XBeach simulation using the Inductiva API, based on a real-world dataset 
that requires significant computational resources.

## Objective
The goal is to demonstrate how to run the `Galveston Island` use case from the [GRIIDC repository](https://data.griidc.org/data/HI.x833.000:0001) - a research data platform from Texas A&M University-Corpus Christi’s Harte Research Institute for Gulf of Mexico Studies.

## Prerequisites
1. Download the dataset titled "XBeach model setup and results for beach and dune enhancement scenarios on Galveston Island, Texas", available [here](https://data.griidc.org/data/HI.x833.000:0001#individual-files).
2. Navigate to **Files >> XBeach_Model_Runs >> Beach_Nourish_Only >> Input** and download *all* the files in this directory.
3. Save the downloaded files in a local folder named `Beach_Nourish_Only`, maintaining this directory structure:

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
To reduce simulation time, update the `params.txt` file with the following changes:
- Add `single_dir = 0` just after the header (required for XBeach v10+).
- Set `tstop` to `34560` to shorten the simulation duration.

## Run Your Simulation
Below is the code required to run the simulation using the Inductiva API.

In this example, we use a `c2d-highcpu-56` cloud machine featuring 56 virtual CPUs (vCPUs) and a 20 GB data disk.

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
	provider="GCP",
    machine_type="c3d-highcpu-56",
    data_disk_gb=20,
    spot=True,)

# Initialize the Simulator
xbeach = inductiva.simulators.XBeach( \
    version="1.24")

# Run simulation 
task = xbeach.run(
    input_dir="Beach_Nourish_Only",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: Setting `spot=True` enables the use of spot machines, which are available at substantial discounts. 
> However, your simulation may be interrupted if the cloud provider reclaims the machine.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 29/06, 18:19:10      8.781 s
	In Queue                  at 29/06, 18:19:19      42.799 s
	Preparing to Compute      at 29/06, 18:20:02      4.636 s
	In Progress               at 29/06, 18:20:06      4598.671 s
		└> 4598.451 s      /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus xbeach params.txt
	Finalizing                at 29/06, 19:36:45      6.671 s
	Success                   at 29/06, 19:36:52      

Data:
	Size of zipped output:    398.16 MB
	Size of unzipped output:  668.30 MB
	Number of output files:   29

Estimated computation cost (US$): 0.40 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 1 hour and 17 minutes.

## Upgrading to Powerful Machines
One of Inductiva’s key advantages is how easily you can scale your simulations to larger, more powerful machines with minimal code changes. Scaling up simply requires updating the `machine_type` parameter when allocating your cloud machine.

You can upgrade to a next-generation cloud machine, increase the number of vCPUs, or do both!

For example, running the simulation on a machine with more vCPUs, such as the `c2d-highcpu-112`, reduces computation time from 1 hour and 17 minutes to approximately **47 minutes**, with a modest cost increase to US$0.48.

It’s that simple! 🚀
