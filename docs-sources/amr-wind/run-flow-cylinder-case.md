# Flow Around a Circular Cylinder 🌀
This tutorial will simulate a flow over a cylinder case using AMR-Wind. This classic fluid dynamics problem reveals the changes in flow behavior depending on the Reynolds number (Re).

We will cover the `ib_cylinder_Re_300` use case from the test files folder of the [AMR-Wind GitHub repository](https://github.com/Exawind/amr-wind/tree/main/test/test_files/ib_cylinder_Re_300).

We will also demonstrate Inductiva’s ability to efficiently scale this use case, starting with a cloud 
machine equivalent to a typical laptop and then scaling up to a more powerful instance.

<img src="_static/Re100.gif" alt="Demo Animation" width="400"/>  <img src="_static/Re10000.gif" alt="Demo Animation" width="400"/>

## Case Description
The flow around a circular cylinder has been extensively studied because it captures fundamental fluid dynamics phenomena. At very low Reynolds numbers (Re < 5), the flow remains steady and symmetric. As Re increases, the flow begins to separate and forms steady recirculation zones. 

When the Reynolds number surpasses approximately 47, the flow transitions to an unsteady regime. Here, periodic vortex shedding develops, creating the characteristic von Kármán vortex street.

At Re = 300, the flow becomes fully three-dimensional and unsteady. This makes it a popular benchmark for validating CFD solvers, immersed boundary methods, and mesh refinement strategies.

## Run the Circular Cylinder Case

### Prerequisites
Download the required files [here](https://storage.googleapis.com/inductiva-api-demo-files/flow-cylinder-case.zip) and place them in a folder called `SimulationFiles`. 

### Case Modifications
To capture the formation of vortices, the original case setup was adjusted as follows:

* Changed the stopping criterion from `max_steps` to `stop_time` to run the simulation for a physical time of 20 seconds.

```diff
- time.stop_time               =   -10.0     # Max (simulated) time to evolve
+ time.stop_time               =   20.0 
- time.max_step                =   20        # Max number of time steps
+ time.max_step                =   -20 
```

* Reduced the plotting frequency to decrease the number of output files.

```diff 
- time.plot_interval            =  10       # Steps between plot files
+ time.plot_interval            =  400      # Reduced output frequency to limit file size
```

* Added vorticity magnitude as a new derived output.

```diff 
- io.derived_outputs = "components(velocity,0,1)" "components(gp,0,1)"
+ io.derived_outputs = "components(velocity,0,1)" "components(gp,0,1)" "mag_vorticity"
```

## Running the Simulation
Below is the code required to run the simulation using the Inductiva API.

In this example, we're using a `c2d-highcpu-16` cloud machine equipped with 16 virtual CPUs (vCPUs), comparable 
in performance to a typical laptop.

```python
"""AMR-Wind example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True)

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind(\
    version="3.4.1")

# Run simulation
task = amr_wind.run(input_dir="/Path/to/SimulationFiles",
    sim_config_filename="ib_cylinder_Re_300.inp",
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
	Waiting for Input         at 12/06, 13:57:05      0.835 s
	In Queue                  at 12/06, 13:57:06      33.458 s
	Preparing to Compute      at 12/06, 13:57:40      2.688 s
	In Progress               at 12/06, 13:57:42      5685.857 s
		└> 5685.675 s      /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus amr_wind ib_cylinder_Re_300.inp
	Finalizing                at 12/06, 15:32:28      2.17 s
	Success                   at 12/06, 15:32:30      

Data:
	Size of zipped output:    288.14 MB
	Size of unzipped output:  852.06 MB
	Number of output files:   2240

Estimated computation cost (US$): 0.14 US$
```

As you can see in the "In Progress" line, the part of the timeline that
represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 1 hour and 35 minutes (5686 seconds).

## Upgrading to Powerful Machines
One of Inductiva’s key advantages is how easily you can scale your simulations to larger, more powerful machines with minimal code changes. Scaling up simply requires updating the `machine_type` parameter when allocating your cloud machine.

You can upgrade to a next-generation cloud machine, increase the number of vCPUs, or do both!

By repeating the simulation on a **c4-highcpu-16** instance, with the same number of vCPUs but two generations newer, 
the runtime is reduced to **73 minutes**, achieving a **1.97× speedup** at a cost of US$0.37.

Alternivately, switching from a cloud machine equivalent to your laptop (**c2d-highcpu-16**) to a machine with more vCPUs (**c2d-highcpu-112**) reduces computation time from **1 hour and 35 minutes** to just 
**37 minutes**, costing US$0.38. This results in a **2.57x faster** simulation!

For more computationally intensive tasks, the benefits of scaling can be even more significant. 🚀

<br>

> To analyze the simulation data programmatically, Python-based tools like **yt** can be used, enabling 
custom visualizations and data extraction. For step-by-step guidance on creating slice plots and animations, 
be sure to check out our [post-processing yt tutorial](https://inductiva.ai/guides/amr-wind/using-yt).