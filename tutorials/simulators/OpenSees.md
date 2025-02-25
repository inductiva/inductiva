In this guide, we will walk you through setting up and running OpenSees simulations
using the Inductiva API.

We will cover an example code to help you get started with simulations.

# OpenSees

[OpenSees (Open System for Earthquake Engineering Simulation)](https://opensees.berkeley.edu/)
is an open-source software framework for simulating the response of structural
and geotechnical systems subjected to earthquakes and other dynamic loads.

Developed by the Pacific Earthquake Engineering Research (PEER) Center, OpenSees
provides a flexible and extensible platform with scripting capabilities through
Tcl and Python. It supports nonlinear analysis, finite element modeling, and
advanced material behavior, making it widely used in academic research and
engineering practice for performance-based seismic design and analysis.

## Supported Versions  
We currently support the following OpenSees versions:  
- **v2.5.0** – Supports Tcl scripting only.  
- **v3.7.1** – Supports both Python and Tcl.

## Example Code

In the following example, we demonstrate how to run an OpenFAST simulation 
using Inductiva's cloud infrastructure. 

```{literalinclude} ../../inductiva/tests/test_simulators/opensees/opensees.py
:language: python
```

The following example demonstrates how to run an OpenSees simulation using MPI and Python.  

To use Tcl instead, simply set the `interface` argument to `tcl`.  

If you want to run a Python simulation without MPI, remove the `n_vcpus` argument
and explicitly set `use_hwthread=False`. Both of these are MPI-related arguments,
and specifying any MPI argument will result in the simulation running with MPI.

## Advanced Tutorial: Running the `SmallMP` Example  

### Objective  

In this tutorial, you'll learn how to use Inductiva's API to run an OpenSees
simulation. We'll validate our simulator implementation while demonstrating its
ability to scale efficiently on more powerful machines.  

### Prerequisites  

The requirements for this tutorial are minimal. Simply download the input files
from [here](https://github.com/OpenSees/OpenSees/tree/master/EXAMPLES/SmallMP).  

Once you have the simulation files, you're ready to scale your simulations to the cloud.  

### Overview

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
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

### Step 1: Running Your Simulation  

#### Choose the Machine for Your Simulation  

We have a wide range of machines available, but since this is a validation
example from the simulator repository and the simulation is relatively short,
we'll select a cost-effective machine with 4 virtual CPUs.  

```python
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-4")
```

#### Initialize Your Simulator  

For this OpenSees simulation, we need to select the appropriate OpenSees class.
This step is straightforward. However, there are a few options to consider.
First, we'll choose the interface; for this simulation, we'll use `tcl`, though
Python is also an option for other simulations. Additionally, we'll select
OpenSees version 3.7.1 for the simulation.  

```python
# Initialize the Simulator
opensees = inductiva.simulators.OpenSees( \
    interface="tcl",
    version="3.7.1")
```

#### Running the Simulation  

At this point, most of the setup is already complete. Now, we simply need to
execute the simulation using our input files.  

```python
# Run the simulation
task = opensees.run( \
    input_dir="/Path/to/SmallMP",
    sim_config_filename="Example.tcl",
    n_vcpus=4,
    on=cloud_machine)
```

In this step, we specify the `sim_config_filename`, which points to the main
`tcl` file that will execute the simulation. We also ensure the simulation
utilizes all available virtual CPUs by setting `n_vcpus=4`.  


#### Analyzing the Results  

Now, it's time to wait for the simulation to complete. As mentioned earlier,
this simulation should run quickly, as its primary goal is to ensure
the simulator is functioning correctly.  

```python
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

As shown here, our simulation took 28.2 seconds in the `In Progress` stage.  

Although this simulation is short, there's still room for improvement in
reducing the processing time.  

Keep reading to discover how easy it is to scale this simulation and cut the
processing time by more than half.  


### Step 2: Scaling Up Your Simulation  

At Inductiva, we aim to make your life easier. Scaling up your simulation is as
simple as changing just two lines of code from the previous script.  

1. Update the `machine_type` to `c2-standard-16`  
2. Increase the number of virtual CPUs to 16  

Your Python script should now look like this:

```python
"""OpenSees Simulation."""
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2-standard-16")

# Initialize the Simulator
opensees = inductiva.simulators.OpenSees( \
    interface="tcl",
    version="3.7.1")

# Run the simulation
task = opensees.run( \
    input_dir="/Path/to/SmallMP",
    sim_config_filename="Example.tcl",
    n_vcpus=16,
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

By changing just these two lines, you're now running your simulation on a much
more powerful machine.  

As a result of this change, here's the updated simulation summary:

```
inductiva tasks info iryyiw0xjufymc9iwch1esabo

Task status: Success

Timeline:
	Waiting for Input         at 25/02, 16:43:55      1.643 s
	In Queue                  at 25/02, 16:43:57      24.61 s
	Preparing to Compute      at 25/02, 16:44:22      3.44 s
	In Progress               at 25/02, 16:44:25      10.219 s
		└> 10.085 s        /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 16 OpenSeesMP Example.tcl
	Finalizing                at 25/02, 16:44:35      1.234 s
	Success                   at 25/02, 16:44:37      

Data:
	Size of zipped output:    23.23 MB
	Size of unzipped output:  59.07 MB
	Number of output files:   186

Estimated computation cost (US$): 0.0011 US$

Go to https://console.inductiva.ai/tasks/iryyiw0xjufymc9iwch1esabo for more details.
```

With just a two-line change, we've reduced the process time from 28.2 seconds to
10.2 seconds. This is the kind of power Inductiva provides to its users.  

Are you in the testing phase, unsure if your simulation has errors or will converge?
No problem! Start with a cost-effective, slower machine. Once you're confident
in your results, seamlessly scale your simulation to full speed with no friction.
