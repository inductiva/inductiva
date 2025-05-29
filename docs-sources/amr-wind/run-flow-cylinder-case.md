# Flow Around a Circular Cylinder 🌀
This tutorial will simulate a flow over a cylinder case using AMR-Wind. This classic fluid dynamics problem reveals the changes in flow behavior depending on the Reynolds number (Re).

We will cover the `ib_cylinder_Re_300` use case from the test files folder of the [AMR-Wind GitHub repository](https://github.com/Exawind/amr-wind/tree/v3.4.0).

<img src="_static/Re100.gif" alt="Demo Animation" width="400"/>  <img src="_static/Re10000.gif" alt="Demo Animation" width="400"/>

## Case Description
The flow around a circular cylinder has been extensively studied because it captures fundamental fluid dynamics phenomena. At very low Reynolds numbers (Re < 5), the flow remains steady and symmetric. As Re increases, the flow begins to separate and forms steady recirculation zones. 

When the Reynolds number surpasses approximately 47, the flow transitions to an unsteady regime. Here, periodic vortex shedding develops, creating the characteristic von Kármán vortex street [1].

At Re = 300, the flow becomes fully three-dimensional and unsteady. This makes it a popular benchmark for validating CFD solvers, immersed boundary methods, and mesh refinement strategies [2].


[1] *[Kármán vortex street – Wikipedia]*(https://en.wikipedia.org/wiki/K%C3%A1rm%C3%A1n_vortex_street)
[2] *[Williamson, C.H.K. (1996). Vortex dynamics in the cylinder wake. Annual Review of Fluid Mechanics.]*(https://doi.org/10.1146/annurev.fl.28.010196.002401)

## Run the Circular Cylinder Case

### Prerequisites
Download the required files [here](https://github.com/Exawind/amr-wind/tree/main/test/test_files/ib_cylinder_Re_300) and place them in a folder called `SimulationFiles`. 

### Case Modifications
To improve the visibility of vortex shedding and optimize computational efficiency, the original case setup was modified with the following changes:

* Changing stopping criteria from `max_steps` to `stop_time`, so that the simulation runs for a physical time of 10s.
```diff
- time.stop_time               =   -10.0     # Max (simulated) time to evolve
+ time.stop_time               =   10.0 
- time.max_step                =   20        # Max number of time steps
+ time.max_step                =   -20 
```

* Increasing time step size to reduce computation time
```diff 
- time.cfl              =   0.45         # CFL factor
+ time.cfl              =   1.0
```

* Decreasing plotting frequency to reduce number of output files 
```diff 
- time.plot_interval            =  10       # Steps between plot files
+ time.plot_interval            =  100      # Reduced output frequency to limit file size
```

* Adding vorticity magnitude as a new derived output
```diff 
- io.derived_outputs = "components(velocity,0,1)" "components(gp,0,1)"
+ io.derived_outputs = "components(velocity,0,1)" "components(gp,0,1)" "mag_vorticity"
```

*  Changing Re from 100 to 1000 to induce faster formation of vortices
```diff 
- transport.viscosity = 1.0e-3   # Set for Re = D*v/mu = 100;
+ transport.viscosity = 1.0e-4   #Adjusted for Re = 1000
```

* Adjusting the flow domain so that vortices are more visible, and  proportional modification is made on the mesh. Also, to reduce computation time, the mesh refinement is reduced from 2 levels to 1. 

```diff 

- amr.n_cell     = 64 64 16   # Grid cells at coarsest AMRlevel
+ amr.n_cell     = 128 64 8 # Doubled x-divisions to match domain size                                                                                         
- amr.max_level = 2
+ amr.max_level = 1

- geometry.prob_lo        =   -0.5 -0.5 -0.125
+ geometry.prob_lo        =   -0.3 -0.5 -0.0625 # Cylinder offset to extend wake region
- geometry.prob_hi        =    0.5  0.5  0.125
+ geometry.prob_hi        =    1.7  0.5  0.0625  
```
* Associated modifications have to be performed on the `static_box.refine` file as well.  Here, we remove the 2nd level of mesh refinement to reduce computation time

```diff

- 2 # Number of levels of refinement
+ 1
  1 # Number of refinement boxes in the first level
 -0.125 -0.125 -0.125 0.5 0.125 0.125 
- 1 # Number of refinement boxes in the second level
- -0.0625 -0.0625 -0.125 0.0625 0.0625 0.125 
```


## Running the Simulation
Here is the code required to run the simulation using the Inductiva API:

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
    version="2.4.0")

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

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

In this example, we're using a cloud machine (`c2d-highcpu-56`) equipped with 56 virtual CPUs.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

<include task print summary>

As you can see in the "In Progress" line, the part of the timeline that
represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 386.5 s.

To analyze the simulation data programmatically, Python-based tools like **yt** can be used, enabling custom visualizations and data extraction. For step-by-step guidance on creating slice plots and animations, be sure to check out our[post-processing yt tutorial](https://inductiva.ai/guides/amr-wind/using-yt).

## Scaling Up the Simulation
One of the key advantages of using Inductiva is the ease with which you can scale your simulations to larger, more powerful machines with minimal changes to your code. Scaling up simply involves updating the `machine_type` parameter when creating your MachineGroup.

Rather than increasing the number of vCPUs on the existing c2d cloud machine, we chose to rerun the simulation on a cloud GPU machine to leverage accelerated hardware. In the following code example, we use the `g2-standard-4` instance type.

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="g2-standard-4",
    data_disk_gb=50,
    spot=True)

# Download the input files into a folder
input_dir = "/Path/to/SimulationFiles"

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind( \
    version="3.4.1")

# Run simulation
task = amr_wind.run(input_dir=input_dir,
    sim_config_filename="ib_cylinder_Re_300.inp",
    n_vcpus=1,
    on=cloud_machine)
```

<table both results and compare them in a sentence>

May your residuals drop fast and your vortices stay coherent. Happy simulating!
