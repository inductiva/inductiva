# Run Your First Simulation
This tutorial will show you how to run XBeach simulations using the Inductiva API. 

We will cover the `Bijleveld_surfbeat_s200` test case from the [official XBeach Subversion repository](https://svn.oss.deltares.nl/repos/xbeach/) to help you get started with simulations.


## Prerequisites
Download all of the required files [here](https://svn.oss.deltares.nl/repos/xbeach/testcases/Wong2016/Bijleveld_surfbeat_s200/) and place them into a folder named `SimulationFiles`. Then, you’ll be ready to send your simulation to the Cloud.

## Adjust Simulation Parameters (Optional)
For the purposes of this tutorial, we'll shorten the simulation time by a factor of 10. To do this, modify the `tstop` parameter in the `params.txt` file to 1800.

## Running an XBeach Simulation
Here is the code required to run an XBeach simulation using the Inductiva API:

```python
"""XBeach example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-standard-30",
	spot=True)

# Initialize the Simulator
XBeach = inductiva.simulators.XBeach( \
    version="1.24")

# Run simulation
task = xbeach.run(input_dir="/Path/to/SimulationFiles",
    sim_config_filename="params.txt",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt the code for this or any other use case, simply replace `input_dir` with the path to your XBeach files before executing it in a Python script.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
    Waiting for Input         at 09/04, 15:17:40      1.219 s
    In Queue                  at 09/04, 15:17:42      45.35 s
    Preparing to Compute      at 09/04, 15:18:27      1.288 s
    In Progress               at 09/04, 15:18:28      230.371 s
        └> 230.219 s       /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus xbeach params.txt
    Finalizing                at 09/04, 15:22:19      5.173 s
    Success                   at 09/04, 15:22:24      

Data:
    Size of zipped output:    131.96 MB
    Size of unzipped output:  175.41 MB
    Number of output files:   12

Estimated computation cost (US$): 0.019 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 230.4 seconds (approximately 4 min).

It's that simple!
