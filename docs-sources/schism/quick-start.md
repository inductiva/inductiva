# Run Your First Simulation
This tutorial will show you how to run SCHISM simulations using the Inductiva API. 

We will cover the `Test_HydraulicStruct` use case from the verification tests collection from the [official SCHISM documentation](https://schism-dev.github.io/schism/master/getting-started/test_suite.html) to help you get started with simulations.


## Prerequisites
Download the required files [here](https://columbia.vims.edu/schism/schism_verification_tests/Test_HydraulicStruct/) and place them in a folder called `SimulationFiles`. Then, you’ll be ready to send your simulation to the Cloud.

## Running an SCHISM Simulation
Here is the code required to run a SCHISM simulation using the Inductiva API:

```python
"""SCHISM example"""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-standard-4",
	spot=True)

# Initialize the Simulator
schism = inductiva.simulators.SCHISM(\
    version="5.11.0")

# Run simulation
task = schism.run(input_dir="/Path/to/SimulationFiles",
    n_vcpus=3,
    num_scribes=2,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other SCHISM simulations, replace `input_dir` with the
path to your SCHISM input files.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 08/04, 14:50:11      1.405 s
	In Queue                  at 08/04, 14:50:12      58.599 s
	Preparing to Compute      at 08/04, 14:51:11      4.892 s
	In Progress               at 08/04, 14:51:15      104.461 s
		├> 1.07 s          mkdir -p outputs
		└> 103.196 s       /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 3 /schism/build/bin/pschism 2
	Finalizing                at 08/04, 14:53:00      16.287 s
	Success                   at 08/04, 14:53:16      

Data:
	Size of zipped output:    540.46 MB
	Size of unzipped output:  705.32 MB
	Number of output files:   31

Estimated computation cost (US$): 0.0013 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 104.5 seconds (approximately 2 minutes).

It's that simple!
