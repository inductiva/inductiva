This guide will walk you through setting up and running OpenSees simulations
using the Inductiva API. 

We will cover two example codes to help you get started with simulations:
- Running OpenSees scripts written in Tcl
- Running OpenSees scripts written in Python (OpeenSeesPy)

# OpenSees

[OpenSees](https://opensees.berkeley.edu/) (Open System for Earthquake Engineering Simulation)
is a software framework designed for the development of sequential, parallel
and grid-enabled finite element applications in earthquake engineering. It
allows users to simulate the response of structural and geotechnical systems
subjected to earthquakes and other hazards using scripts written in either Tcl
or Python.

The software provides advanced capabilities for modelling and analyzing the
non-linear response of systems, offering a wide range of material models,
elements and solution algorithms.


## Supported Versions
We currently support the following OpenSees versions:
- **v2.5.0** - Supports Tcl scripting only.
- **v3.7.1** - Supports Python and Tcl scritping.

## Running OpenSees scripts written in Tcl

### Objective

This tutorial will show you how to run an OpenSees simulation using its Tcl
interface, using the `SmallMP` use case available in the
[OpenSees GitHub repository](https://github.com/OpenSees/OpenSees).

We will also demonstrate Inductiva's ability to efficiently scale this use case
on a more powerful machine.

### Prerequisites  

The requirements for this tutorial are minimal. Simply download the input files
from [here](https://github.com/OpenSees/OpenSees/tree/master/EXAMPLES/SmallMP).  

Once you have the simulation files, you're ready to scale your simulations to
the Cloud.

### Running Your Simulation

Here is the code required to run an OpenSees simulation using the Inductiva API.

```python
"""OpenSees Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-4")

# Initialize the Simulator
opensees = inductiva.simulators.OpenSees( \
    interface="tcl",
    version="3.7.1")

# Run simulation
task = opensees.run( \
    input_dir="/Path/to/SmallMP",
    sim_config_filename="Example.tcl",
    n_vcpus=4,
    use_hwthread=True,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

To adapt it for this or any other use case, simply replace `input_dir` with the
path to your OpenSees files and specify the `sim_config_filename` before running
it in a Python script.

Once the simulation is complete, we terminate the machine, download the results
and print a summary of the simulation as shown below.

```
inductiva tasks info gafdcf5t0zpkft4sxmubo30q1

Task status: Success

Timeline:
	Waiting for Input         at 25/02, 16:35:31      1.995 s
	In Queue                  at 25/02, 16:35:33      20.562 s
	Preparing to Compute      at 25/02, 16:35:54      3.996 s
	In Progress               at 25/02, 16:35:58      28.221 s
		└> 28.09 s         /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 4 OpenSeesMP Example.tcl
	Finalizing                at 25/02, 16:36:26      1.339 s
	Success                   at 25/02, 16:36:27      

Data:
	Size of zipped output:    23.23 MB
	Size of unzipped output:  59.07 MB
	Number of output files:   186

Estimated computation cost (US$): 0.00066 US$

Go to https://console.inductiva.ai/tasks/gafdcf5t0zpkft4sxmubo30q1 for more details.
```

The core computation time of our simulation was 28.2 seconds, as can be seen in
the line `In Progress at 25/02, 16:35:58 28.221 s`. This part of the timeline
represents the actual execution of the simulation.

Although it's short, there's still room for improvement to reduce the processing
time.

### Scaling Up Your Simulation  

Scaling up your simulation is as simple as changing just two lines of code:

1. Update the `machine_type` to a 16 vCPU machine (`c2-standard-16`)
2. Set the number of `n_vcpus` to 16  

By increasing the number of vCPUs, we've reduced the processing time from 28.2 
to 10.2 seconds.

Here are the results of running the same simulation on a few machines:

|  Machine Type  | Virtual CPUs |     Time     | Estimated Cost |
|:--------------:|:------------:|:------------:|:--------------:|
|  c2-standard-4 |       4      | 28.2 seconds | 0.00066 US$    |
|  c2-standard-8 |       8      | 15.2 seconds | 0.00094 US$    |
| c2-standard-16 |      16      | 10.2 seconds | 0.0011 US$     |

Still in the testing phase? No problem! Just skip this step for now and start
with a machine with fewer vCPUs. Once you're satisfied with your results, you
can seamlessly scale your OpenSees simulation.

## Running OpenSees scripts written in Python (OpeenSeesPy)

This tutorial will show you how to run an OpenSees simulation using its Python
interface (OpenSeesPy), using one of the Python examples referenced in
[OpenSees GitHub repository](https://github.com/OpenSees/OpenSees).

For the purposes of this tutorial, we decided to download 
[this input file](https://github.com/OpenSees/OpenSees/blob/master/EXAMPLES/ExamplePython/example_mpi_paralleltruss_explicit.py).

Once you have downloaded the file, you are ready to execute the code required to
run an OpenSees simulation.

```python
"""OpenSees Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-4")

# Initialize the Simulator
opensees = inductiva.simulators.OpenSees( \
    interface="python",
    version="3.7.1")

# Run simulation
task = opensees.run( \
    input_dir="/Path/to/ExamplePython",
    sim_config_filename="example_mpi_paralleltruss_explicit.py",
    n_vcpus=4,
    use_hwthread=True,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

All you need to do is change the `interface` parameter to `python` to match the
file type of your OpenSeesPy use case. It's that simple!
