# Flow around a circular cylinder 🌀
This tutorial will show you how to run AMR-Wind simulations using the Inductiva API. 

This tutorial will simulate a flow over a cylinder case using amr-wind. The `ib_cylinder_Re_300` use case from the test files folder of the [AMR-Wind GitHub repository](https://github.com/Exawind/amr-wind/tree/) to help you get started with simulations.


  <img src="_static/Re100.gif" alt="Demo Animation" width="400"/>  <img src="_static/Re10000.gif" alt="Demo Animation" width="400"/>

## Case description
The flow around a circular cylinder is a classic problem in fluid dynamics, with behavior that changes dramatically depending on the Reynolds number (Re), which compares inertial and viscous forces. At very low Re (< 5), the flow remains steady and symmetric. As Re increases, flow separation occurs, followed by the formation of steady recirculation zones . Beyond Re ≈ 47, the flow becomes unsteady, exhibiting periodic vortex shedding that leads to the formation of a von Kármán vortex street [1]. At Re = 300, the flow is fully three-dimensional and unsteady, making it a widely used benchmark for validating CFD solvers, immersed boundary methods, and mesh refinement strategies [2].

### References
[1] [Kármán vortex street – Wikipedia](https://en.wikipedia.org/wiki/K%C3%A1rm%C3%A1n_vortex_street)

[2] [Williamson, C.H.K. (1996). Vortex dynamics in the cylinder wake. Annual Review of Fluid Mechanics.](https://doi.org/10.1146/annurev.fl.28.010196.002401)

## Prerequisites
Download the required files [here](https://github.com/Exawind/amr-wind/tree/main/test/test_files/ib_cylinder_Re_300) and place them in a folder called `SimulationFiles`. Then, you’ll be ready to send your simulation to the Cloud.

## Case Modifications

To improve the visibility of vortex shedding and optimize computational efficiency, the original case setup was modified with the following changes:

*    Increased domain length in the x-direction to 2.0 units for better capture of the wake development
*    Simplified mesh refinement, including removal of the second refinement level and corresponding updates in `static_box.refine`
*    Adjusted simulation time and output frequency, increasing total simulation time while reducing save frequency
*    Updated flow parameters by modifying viscosity to achieve a Reynolds number of 1000
*    Added `mag_vorticity` to `io.derived_outputs` to visualize vortex structures

```diff
## ib_cylinder_Re_300.inp ##

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
- time.stop_time               =   -10.0     # Max (simulated) time to evolve
+ time.stop_time               =   10.0      # Switched stop condition to physical time
- time.max_step                =   20        # Max number of time steps
+ time.max_step                =   -20       # Making the term negative to deactivate it

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   -0.05        # Use this constant dt if > 0
- time.cfl              =   0.45         # CFL factor
+ time.cfl              =   1.0

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
- time.plot_interval            =  10       # Steps between plot files
+ time.plot_interval            =  100      # Reduced output frequency to limit file size

time.checkpoint_interval      =  -1       # Steps between checkpoint files

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
ConstValue.density.value = 1.0
ConstValue.velocity.value = 1.0 0.0 0.0

io.output_default_variables = 0
io.outputs = density p
- io.derived_outputs = "components(velocity,0,1)" "components(gp,0,1)"
+ io.derived_outputs = "components(velocity,0,1)" "components(gp,0,1)" "mag_vorticity"
incflo.use_godunov = 1
incflo.diffusion_type = 2
incflo.godunov_type = "weno_z"
incflo.do_initial_proj = 1
incflo.initial_iterations = 3
- transport.viscosity = 1.0e-3   # Set for Re = D*v/mu = 100;
+ transport.viscosity = 1.0e-4   #Adjusted for Re = 1000
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = Laminar

incflo.physics = FreeStream IB
IB.labels = IB1  
IB.IB1.type = Cylinder 
IB.IB1.center = 0.0 0.0 0.0
IB.IB1.radius = 0.05 
IB.IB1.height = 0.25

- amr.n_cell     = 64 64 16   # Grid cells at coarsest AMRlevel
+ amr.n_cell     = 128 64 8 # Doubled x-divisions to match domain size
tagging.labels = sr                                                                                                                
tagging.sr.type = CartBoxRefinement                                                                                                                
tagging.sr.static_refinement_def = static_box.refine                                                                                             
- amr.max_level = 2
+ amr.max_level = 1

- geometry.prob_lo        =   -0.5 -0.5 -0.125
+ geometry.prob_lo        =   -0.3 -0.5 -0.0625 # Cylinder offset to extend wake region
- geometry.prob_hi        =    0.5  0.5  0.125
+ geometry.prob_hi        =    1.7  0.5  0.0625  
geometry.is_periodic    =   0   0   1   # Periodicity x y z (0/1)

# Boundary conditions
xlo.type = "mass_inflow"
xlo.density = 1.0
xlo.velocity = 1.0 0.0 0.0
xhi.type = "pressure_outflow"
ylo.type =   "slip_wall"
yhi.type =   "slip_wall"

incflo.verbose          =   0          # incflo_level
nodal_proj.verbose = 0

nodal_proj.mg_rtol = 1.0e-12
nodal_proj.mg_atol = 1.0e-12
mac_proj.mg_rtol = 1.0e-12
mac_proj.mg_atol = 1.0e-12

```

```diff
## static_box.refine ##

- 2 # Number of levels of refinement
+ 1
  1 # Number of refinement boxes in the first level
 -0.125 -0.125 -0.125 0.5 0.125 0.125 #Defining the boundaries of refinement box <xlo ylo zlo xhi yhi zhi>
- 1 # Number of refinement boxes in the second level
- -0.0625 -0.0625 -0.125 0.0625 0.0625 0.125 #Defining the boundaries of refinement box <xlo ylo zlo xhi yhi zhi>
```

## Running an AMR-Wind Simulation
Here is the code required to run a AMR-Wind simulation using the Inductiva API:

```python
"""AMR-Wind example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-56",
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

When the simulation is complete the machine group will be automatically terminated by the `terminate()` command. 

Once the simulation is done, the output files will automotically downloaded to the local machine (if the simulation was run on terminal). If not the files can be downloaded from the 'tasks' tab in the console. 

## Post-Processing

After running the simulation, AMR-Wind generates output files for each time step, typically of the format `plt00000`. These plotfiles contain the simulation data and can be visualized directly using tools like [ParaView](https://www.paraview.org/), which supports AMReX plotfiles natively. To analyze the data programmatically, Python-based tools such as [yt](https://yt-project.org/) can be employed, allowing for custom visualizations and data extraction. Additionally, if the sampling feature is enabled in your simulation setup, selected data can be exported in standard formats like NetCDF, typically saved in the `post_processing/` directory for easy access and compatibility with other tools.

Checkout how you can use `yt` to generate custom slice plots and animations [here](using-yt.md).

## Speedup

The above simulation uses a [C2D machine](https://cloud.google.com/compute/docs/compute-optimized-machines#c2d_series) with 56 CPU cores from Google Cloud. This simulation can be run much faster by increasing the number of CPUs or by running the simulation on a GPU machine. The definition for a GPU machine is provided below.

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
If you're wondering how much faster the simulation can run on different machines—or what the cost trade-offs look like— we have benchmarked a [20-million-cell CFD simulation](input_fine_Re10000.inp) across a variety of CPU and GPU instances. The results provided below offer a practical comparison of run time and cost efficiency across platforms. Whether you're optimizing for speed, budget, or a balance of both, this data can help you decide what kind of machine is best suited for your needs.

## Benchmark Performance Comparison: CPU vs. GPU Configurations 📈

This page summarizes the results of a comprehensive benchmark test comparing execution time, cost, and performance across multiple CPU and GPU configurations. 
The goal is to provide insights into the optimal hardware setup for this workload, balancing speed, resource utilization, and cost-effectiveness. A 20M-cell CFD simulation 
of flow over a cylinder at Re=10,000 was used for this purpose. The input file used for this simulation is provided [here](input_fine_Re10000.inp) for reproducibility.


<table>
    <tr>
        <td>Platform</td>
        <td>Machine Type</td>
        <td>Machines</td>
        <td>Cores</td>
        <td>GPUS</td>
        <td>Duration</td>
        <td>Cost (USD)</td>
    </tr>
    <tr>
        <td rowspan="4">AMD Epyc Milan</td>
        <td>c2d-highcpu-16</td>
        <td>1</td>
        <td>16</td>
        <td>0</td>
        <td>2 hours, 8 minutes</td>
        <td>0.23 US$</td>
    </tr>
    <tr>
        <td>c2d-highcpu-32</td>
        <td>1</td>
        <td>32</td>
        <td>0</td>
        <td>1 hour, 14 minutes</td>
        <td>0.26 US$</td>
    </tr>
    <tr>
        <td>c2d-highcpu-56</td>
        <td>1</td>
        <td>56</td>
        <td>0</td>
        <td>53 minutes, 29 seconds</td>
        <td>0.32 US$</td>
    </tr>
    <tr>
        <td>c2d-highcpu-112</td>
        <td>1</td>
        <td>112</td>
        <td>0</td>
        <td>27 minutes, 48 seconds</td>
        <td>0.33 US$</td>
    </tr>
    <tr>
        <td rowspan="4">Nvidia L4 Cuda Cores 7680;  Tensor Cores 240</td>
        <td>g2-standard-4</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <td>55 minutes, 26 seconds</td>
        <td>0.26 US$</td>
    </tr>
    <tr>
        <td>g2-standard-24</td>
        <td>1</td>
        <td>2</td>
        <td>2</td>
        <td>40 minutes, 43 seconds</td>
        <td>0.54 US$</td>
    </tr>
    <tr>
        <td>g2-standard-48</td>
        <td>1</td>
        <td>4</td>
        <td>4</td>
        <td>25 minutes, 7 seconds</td>
        <td>0.66 US$</td>
    </tr>
    <tr>
        <td>g2-standard-96</td>
        <td>1</td>
        <td>8</td>
        <td>8</td>
        <td>18 minutes, 14 seconds</td>
        <td>0.96 US$</td>
    </tr>
<tr>
        <td rowspan="4">Nvidia a100 Cuda Cores 6912;  Tensor Cores 432</td>
        <td>a2-highgpu-1g</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <td>20 minutes, 42 seconds</td>
        <td>0.51 US$</td>
    </tr>
    <tr>
        <td>a2-highgpu-2g</td>
        <td>1</td>
        <td>2</td>
        <td>2</td>
        <td>24 minutes, 50 seconds</td>
        <td>1.23 US$</td>
    </tr>
    <tr>
        <td>a2-highgpu-4g</td>
        <td>1</td>
        <td>4</td>
        <td>4</td>
        <td>17 minutes, 56 seconds</td>
        <td>1.78 US$</td>
    </tr>
    <tr>
        <td>a2-highgpu-8g</td>
        <td>1</td>
        <td>8</td>
        <td>8</td>
        <td>16 minutes, 45 seconds</td>
        <td>3.33 US$</td>
    </tr>
    <tr>
        <td rowspan="2">Nvidia H100 Cuda Cores 14592</td>
        <td>a3-highgpu-1g</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <td>13 minutes, 59 seconds</td>
        <td>0.58 US$</td>
    </tr>
    <tr>
        <td>a3-highgpu-2g</td>
        <td>1</td>
        <td>2</td>
        <td>2</td>
        <td>16 minutes, 15 seconds</td>
        <td>1.36 US$</td>
    </tr>
</table>

May your residuals drop fast and your vortices stay coherent — happy simulating!