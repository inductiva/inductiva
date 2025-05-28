# Run a Temporal Boundary Layer with Stable Stratification Case
---
*This tutorial was written by* **[Pedro Simões]**(P.SimoesCosta@tudelft.nl) *in collaboration with the* **Inductiva Team**

The numerical simulation of a temporally evolving, stably stratified boundary
layer offers a clear computational sandbox for exploring fundamental fluid
dynamics relevant to atmospheric wind flows.

Understanding boundary layer physics is critical in atmospheric modeling because
wind characteristics, turbulence structure, and mixing processes strongly depend
on boundary layer dynamics. The temporal boundary layer (TBL) setup is useful
for understanding boundary layer turbulence at high Reynolds numbers typical of
atmospheric flows. The higher the Reynolds number of a wind flow boundary layer,
the more "locally parallel" the flow is. The TBL flow consists of a moving bottom
wall with constant velocity, with turbulence being entrained to the flow vertically
and parallel to the wall. This asymptotically approximates the development of
windflow over a surface as it approaches a parallel state, allowing simple
interpretation of turbulent processes.

Stability plays a crucial role in altering boundary layer dynamics. To better understand its effects, let's examine the types of stratification conditions:

- **Neutral stratification** (no vertical temperature gradient) provides a baseline scenario where temperature doesn't affect the wind turbulence dynamics.

<img src="_static/tempField_neutralTDBL_Re1000-13863.gif" alt="Demo Animation"/>

- In **stable stratification**, temperature decreases with the height (that is, denser air near the surface and lighter air above), and buoyancy suppresses vertical movements, reducing turbulent mixing, and overall boundary layer growth.

<img src="_static/tempField_neutralTDBL_Re1000-13863.gif" alt="Demo Animation"/>

- Conversely, in **unstable stratification** (cooler air above warmer air), buoyancy amplifies vertical motions, enhancing turbulence and thickening the boundary layer.

Simulating these different stratification conditions provides insight into how buoyancy and turbulence interact to shape the wind boundary layer.

This tutorial focuses on simulating **stable stratification**, where warmer air overlies cooler air, as it closely reflects common atmospheric conditions - particularly at night due to radiative cooling of the land surface or during stable weather patterns.

## Simulate Stable Stratification

### Requirements
To run this use case, begin by creating the `input.nml` with the
following contents:

```
&dns
ng(1:3) = 768, 768, 512
l(1:3) = 150, 75, 80
gtype = 2, gr = 2.
cfl = 0.95, dtmax = 1.e5
visci = 1000.
inivel = 'tbl'
is_wallturb = F
nstep = 40000, time_max = 2000., tw_max = 0.1
stop_type(1:3) = T, F, F
restart = F, is_overwrite_save = T, nsaves_max = 0
icheck = 10, iout0d = 100, iout1d = 100, iout2d = 200, iout3d = 20000, isave = 5000
cbcvel(0:1,1:3,1) = 'P','P',  'P','P',  'D','N'
cbcvel(0:1,1:3,2) = 'P','P',  'P','P',  'D','N'
cbcvel(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
bcvel(0:1,1:3,1) =  0.,0.,   0.,0.,   1.,0.
bcvel(0:1,1:3,2) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,3) =  0.,0.,   0.,0.,   0.,0.
bcpre(0:1,1:3)   =  0.,0.,   0.,0.,   0.,0.
bforce(1:3) = 0., 0., 0.
is_forced(1:3) = F, F, F
velf(1:3) = 0., 0., 0.
gacc(1:3) = 0., 0., -1.
nscal = 1
dims(1:2) = 0, 0, ipencil_axis = 3
/

&scalar
iniscal(:)             = 'tbl'
alphai(:)              = 710.  ! = Pr/nu = 500 * 0.71
beta                   = 0. ! = Gr/Re_D^2, Gr = 1e3, Re_D=500
cbcscal(0:1,1:3,:)     = 'P'  ,'P' ,  'P','P',  'D','N'
bcscal(0:1,1:3,:)      =  0.,0. ,   0.,0. ,   1.,0.
is_sforced(:)          = F
scalf(:)               = 0.
is_boussinesq_buoyancy = T
/

&cudecomp
cudecomp_t_comm_backend = 0, cudecomp_is_t_enable_nccl = T, cudecomp_is_t_enable_nvshmem = T
cudecomp_h_comm_backend = 0, cudecomp_is_h_enable_nccl = T, cudecomp_is_h_enable_nvshmem = T
/

&numerics
is_impdiff = F, is_impdiff_1d = T
is_poisson_pcr_tdma = F
/

&other_options
is_debug = T, is_timing = T
/
```

To shorten the simulation time, we will run only 1% of the original simulation.
To do this, update the following line:

```
nstep = 40000, time_max = 2000., tw_max = 0.1  
```

to:

```
nstep = 400, time_max = 20., tw_max = 0.1  
```

This adjustment significantly reduces the simulation duration while maintaining the overall configuration. If desired, the full simulation can be run by keeping the original values.

### Running the Simulation
Here is the code required to run the simulation using the Inductiva API:

```python
"""CaNS example."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
        provider="GCP",
        machine_type="a3-highgpu-1g",
        zone="europe-west4-b",
        data_disk_gb=100,
        spot=True)

# Initialize the Simulator
cans = inductiva.simulators.CaNS(\
    version="3.0.0")

# Run simulation
task = cans.run(input_dir="/Path/to/SimulationFiles",
    sim_config_filename="input.nml",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()

```

> **Note**: `spot` machines are a lot cheaper but may be terminated by the provider if necessary.

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 22/05, 09:13:45      0.818 s
	In Queue                  at 22/05, 09:13:46      57.083 s
	Preparing to Compute      at 22/05, 09:14:43      22.052 s
	In Progress               at 22/05, 09:15:05      493.001 s
		├> 2.083 s         mkdir -p data
		└> 490.643 s       /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus --np 1 cans input.nml
	Finalizing                at 22/05, 09:23:18      404.624 s
	Success                   at 22/05, 09:30:02      

Data:
	Size of zipped output:    15.59 GB
	Size of unzipped output:  19.37 GB
	Number of output files:   42

Estimated computation cost (US$): 0.65 US$
```

As you can see in the "In Progress" line, the part of the timeline that
represents the actual execution of the simulation, 
the core computation time of this simulation was approximately 8 minutes and 13 seconds.

It's that simple!

### Scalling Up the Simulation
One of the benefits of using Inductiva is the ability to scale simulations to bigger and faster machines with minimal code changes. In this case, only the `machine_type` argument needs to be updated during MachineGroup creation.

Below are the results of the same simulation executed on multiple machines:

| Machine Type    | vCPUs | GPU                | GPU Count | Duration          | Cost      | Speed-up |
| --------------- | ----- | ------------------ | --------- | ----------------- | --------- | -------- |
| c3d-highcpu-90  | 90    | -                  | 0         | 28 minutes 11 sec | 0.49 US$ | 0.28x    |
| c3d-highcpu-180 | 180   | -                  | 0         | 17 minutes 16 sec | 0.64 US$ | 0.46x    |
| c3d-highcpu-360 | 360   | -                  | 0         | 10 minutes 20 sec | 0.89 US$ | 0.78x    |
| a3-highgpu-1    | 26    | NVIDIA H100 (80GB) | 1         | 8 minutes 1 sec   | 0.63 US$ | 1.00x    |
| a3-highgpu-2    | 52    | NVIDIA H100 (80GB) | 2         | 9 minutes 1 sec   | 1.25 US$ | 0.89x    |
| a3-highgpu-4    | 104   | NVIDIA H100 (80GB) | 4         | 6 minutes 5 sec   | 2.05 US$ | 1.32x    |
| a3-highgpu-8    | 208   | NVIDIA H100 (80GB) | 8         | 6 minutes 8 sec   | 3.99 US$ | 1.30x    |


It's important to note that we only ran 1% of the original CaNS simulation.
Any performance improvements observed here would be even more pronounced when
running the full simulation. 

These results also highlight that scaling
up does not always lead to linear speedups. In some cases, it can even lead to
longer runtimes. This is often due to overheads such as communication between
GPUs or underutilization of available resources.

Careful benchmarking is essential to find the right balance between speed,
efficiency, and cost.

Stay tunned for more!