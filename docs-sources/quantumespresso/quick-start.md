# Run Your First Simulation
This tutorial will show you how to run Quantum ESPRESSO simulations using the Inductiva API. 

To help you get started with simulations, we will explore a use case from this [Atomsk tutorial](https://atomsk.univ-lille.fr/tutorial_QE.php), which guides us in creating a unit cell for fcc aluminium.

## Prerequisites
Copy and create your `Al.pw` file, placing it in a designated folder. Then, make the following adjustments:
- Remove the line of code that specifies the `pseudo` directory path: `pseudo_dir = '/home/user/espresso/pseudo/'`.
- Under the section labeled *ATOMIC_SPECIES*, update the pseudopotential file: Change `"Al 26.982 Al.fixme.upf"` to `"Al 26.982 Al.pbe-nl-rrkjus_psl.1.0.0.UPF`.

You're ready to send your simulation to the Cloud!

## Running a Quantum ESPRESSO Simulation
Here is the code required to run a Quantum ESPRESSO simulation using the Inductiva API:

```python
"""Quantum ESPRESSO example."""
import inductiva
from inductiva.commands import MPIConfig, Command

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-standard-4",
    spot=True)

mpi_config = MPIConfig( \
    version="4.1.6",
    np=4,
    use_hwthread_cpus=True)

# List of commands to run
commands = [
    Command("pw.x -i Al.pw", mpi_config=mpi_config),
]

# Initialize the Simulator
qe = inductiva.simulators.QuantumEspresso(\
    version="7.4.1")

# Run simulation
task = qe.run( \
    input_dir="/Path/to/SimulationFiles",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

To adapt this script for other Quantum ESPRESSO simulations, replace `input_dir` with the
path to your Quantum ESPRESSO input files.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 03/04, 16:48:46      0.846 s
	In Queue                  at 03/04, 16:48:46      16.256 s
	Preparing to Compute      at 03/04, 16:49:03      8.963 s
	In Progress               at 03/04, 16:49:12      4.162 s
		└> 4.066 s         /opt/openmpi/4.1.6/bin/mpirun --np 4 --use-hwthread-cpus pw.x -i Al.pw
	Finalizing                at 03/04, 16:49:16      0.459 s
	Success                   at 03/04, 16:49:16      

Data:
	Size of zipped output:    607.28 KB
	Size of unzipped output:  1.48 MB
	Number of output files:   11

Estimated computation cost (US$): 0.00015 US$
```

As you can see in the "In Progress" line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was 4.1 seconds.

It's that simple!
