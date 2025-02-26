This guide will walk you through setting up and running OpenSees simulations
using the Inductiva API. 

We will cover two example codes to help you get started with simulations:
- Running OpenSees scripts written in Tcl
- Running OpenSees scripts written in Python (OpeenSeesPy)

# OpenSees

[OpenSees (Open System for Earthquake Engineering Simulation)](https://opensees.berkeley.edu/)
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
interface, using the `SmallMP` use case present in the
[OpenSees GitHub repository](https://github.com/OpenSees/OpenSees).

We will also demonstrate Inductiva's ability to efficiently scale this use case
on a more powerful machine.

### Prerequisites  

The requirements for this tutorial are minimal. Simply download the input files
from [here](https://github.com/OpenSees/OpenSees/tree/master/EXAMPLES/SmallMP).  

Once you have the simulation files, you're ready to scale your simulations to
the cloud.  

### Running Your Simulation

Here's the code you'll be working on as we progress through the tutorial.
Don't worry if it doesn't all make sense right now; everything will become
clearer in the upcoming steps.

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

Once the simulation finishes, we terminate the machine, download the results,
and print a summary of the simulation.

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

As shown here, the core computing time of our simulation took 28.2 seconds, as
seen at the `In Progress   at 25/02, 16:35:58  28.221 s` line. This part of the
timeline represents the actual running of the simulation.

Although this simulation is short, there's still room for improvement in
reducing the processing time.

### Scaling Up Your Simulation  

In order to scale your simulation you just need to change a couple of lines on
your original script.

1. Update the `machine_type` to `c2-standard-16`  
2. Increase the number of `n_vcpus` to 16  

By changing just these two lines, you're now running your simulation on a much
more powerful machine.

With just a two-line change, we've reduced the process time from 28.2 seconds to
10.2 seconds.

Here are the results of running the same simulation on a few machines:

|  Machine Type  | Virtual CPUs |     Time     | Estimated Cost |
|:--------------:|:------------:|:------------:|:--------------:|
|  c2-standard-4 |       4      | 28.2 seconds | 0.00066 US$    |
|  c2-standard-8 |       8      | 15.2 seconds | 0.00094 US$    |
| c2-standard-16 |      16      | 10.2 seconds | 0.0011 US$     |

Are you in the testing phase, unsure if your simulation has errors or will converge?
No problem! Start with a cost-effective, slower machine. Once you're confident
in your results, seamlessly scale your simulation to full speed with no friction.

## Running OpenSees scripts written in Python (OpeenSeesPy)

In this example, we demonstrate how to run an OpenSees simulation using its
Python interface (OpenSeesPy). This tutorial shares many similarities with the
previous one, so we'll skip over some parts and focus on the key differences.

Before running the simulation, you'll need to download the `ExamplePython` from
the [OpenSees official repository](https://github.com/OpenSees/OpenSees/tree/master/EXAMPLES/ExamplePython).

With the simulation files downloaded you can run your simulation with the
following Python script.

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

### What's Changed

The code in this example is largely the same as the previous one, but with two
key differences:

1. When initializing the simulator, we now specify `interface="python"`.
2. The `sim_config_filename` now points to a Python file: `example_mpi_paralleltruss_explicit.py`.

That's it! All you have to do is change the `interface` parameter to match the
file type of your OpenSees use case.
