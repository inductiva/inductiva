# Reef3D

REEF3D is an open-source hydrodynamics framework with a focus on coastal, marine 
and hydraulic engineering flows. Tailor-made multiphysics solvers are available 
for a range of relevant problems (e.g. sediment transport or floating body
dynamics). The modular programming approach allows the framework to incorporate
a range of different flow solvers which together represent all relevant length
scales. Depending on the wave or flow conditions, the following optimized
hydrodynamic modules are available:

- **REEF3D::CFD** solves the Navier-Stokes equations in three dimensions. For 
near-field simulations with a complex free surface pattern,  it uses a two-phase 
flow approach with the level set method for interface capturing.
- **REEF3D::FNPF** is a three-dimensional fully nonlinear potential flow solver. 
It is massively parallelized and can be used to create large-scale
phase-resolved sea states at all water depths.
- **REEF3D::SFLOW** is a depth-averaged model, solving the non-hydrostatic
shallow water equations ideal for near-shore hydrodynamics and river flow.

## Running a simulation

Reef3D in Inductiva API executes two sequential steps:

- the meshing with **DiveMESH**;
- the simulation with **Reef3D**. 

Each step is configured with input files, `control.txt` and `ctrl.txt`,
respectively. Other files may be used to inform the simulator about the grid,
geographical data or wave information. Reef3D has strict naming policies for
each file and we recommend users to follow
[their guidelines](https://reef3d.wordpress.com/user-guide/). 

To run the simulation users pass a folder with the above files and all others 
required to run the simulation. The folder is then uploaded to Inductiva API and 
the simulation is executed. The output files are then downloaded to the user's 
machine. 

The parallelization of the simulation is handled automatically by Inductiva API, 
based on the number of cores available in the machine.

**General Arguments:**
- `on`: set the machines where the simulations will run. Check
[here](https://tutorials.inductiva.ai/intro_to_api/computational-infrastructure.html#available-computational-resources) 
for further detail. If not selected the simulations will be picked-up by a
default pool shared by everyone.
- `storage_dir`: set the directory where the output files will be stored in the 
cloud. If not selected the output files will be stored in a folder named with
the  task id of the simulation.
- `n_vcpus`: number of virtual CPUs / threads that will be used to configure the
MPI parallism. This number needs to be set consistently with parameter
```M 10``` to be set in both `control.txt` and `ctrl.txt` configurations files.


For further information on handling the task of the simulation see
[here](https://tutorials.inductiva.ai/intro_to_api/tasks.html).

## Example

Here, we follow the tutorial with
[regular wave propagation](https://github.com/REEF3D/REEF3D/tree/ed0c8d7a6110892706357f72e0404bd63034efa5/Tutorials/REEF3D_FNPF/9_1%20Regular%20Wave%20Propagation)
from Reef3D repository.

```python
import inductiva

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip", unzip=True
)

reef3d = inductiva.simulators.REEF3D()

task = reef3d.run(input_dir=input_dir)

task.wait()
task.download_outputs()
```

## A slighly more advanced example
Let's now run a more advanced example, one that will also require a lot more
compute power, and will illustrate more advanced features of the API. More
specifically, we will use the API to run the "3D Dam Break Scenarion with
Obstacle" that can be found in 
[Reef3D tutorials](https://github.com/REEF3D/REEF3D/tree/master/Tutorials/REEF3D_CFD/10_2%203D%20Dam%20Break%20with%20Obstacle). As explained before, there are
two files that are required to configure the simulation:

* [control.txt](https://github.com/REEF3D/REEF3D/blob/master/Tutorials/REEF3D_CFD/10_2%203D%20Dam%20Break%20with%20Obstacle/control.txt)
* [ctrl.txt](https://github.com/REEF3D/REEF3D/blob/master/Tutorials/REEF3D_CFD/10_2%203D%20Dam%20Break%20with%20Obstacle/ctrl.txt)

Let's start by downloading both of these files and adding them to a folder named
"10_2_3D_Dam_Break_with_Obstacle" inside our working folder. Just download each
file directly from GitHub (use "Download raw file" option) and move them to
the "10_2_3D_Dam_Break_with_Obstacle" local folder (that you also need to create
locally).

Before we proceed, let's inspect the files to check three Reed3D parameters that
are important to understand before we configure our simulation run. These
parameters are: ```N 41```, ```P 30 ``` and ```M 10```:

### control.txt (for DiveMESH)
```
C 11 21
C 12 21
C 13 21
C 14 21
C 15 21
C 16 21

B 1 0.025
B 10 0.0 2.0 0.0 1.0 0.0 1.0
O 10 1.2 1.4 0.4 0.6 0.0 1.0

M 10 4    <---- defines the nr. of processors for parallel computations (4)
```

### ctrl.txt (for Reef3D)
```
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 54 0.5
F 56 0.7
N 40 3
N 41 25.0    <---- set the maximum modeled time (25 seconds).
N 45 50000
N 47 0.2
M 10 4    <---- defines the nr. of processors for parallel computations (4)
P 10 1
P 30 0.01    <---- defines the rate of paraview results (1 frame per 0.01 s)
T 10 0
W 22 -9.81
```

Observe that parameter ```M 10```, which controls the degree of parallism, is
set to 4 threads/vCPUs, a very low number. Depending on the number of vCPUs we
effectively wish to use for running the simulation, we will need to manually
change  ```M 10``` on **both files** to match the specs of corresponding the VM. 
Because this is already a reasoably heavy simulation, we will be using GCP VMs
of the c3d family, supported by last-generation AMD chips. More specifically,
we will be using ```c3d-standard-90``` machines with 90 vCPU. So, we will be
setting ```M 10 90``` on both files.

Also, Reef3D produces a huge amount of data. As it is currently configured, this
simulation would produce several dozen gigabytes of data. To reduce that amount
of data produced, we can reduce the maximum modeled time (```N 41```) and
the rate at which Paraview data is being produced (```P 10```).

So, for a reducing the amount of data generated, we will be reducing the rate
at which Paraview information is going to be generated to "onyy" 25 frames per
second, that is, we will set ```P 10 0.04```. 

We are not going to reduce the maximum modeled time but, instead, we will add
extra disk capacity to ensure we have enough space to write all the data will be
generated. Therefore, we will request our machine to be equipped
with **20GB** just for data, using the ```data_disk_gb``` parameter of the 
```inductiva.resources.MachineGroup``` class.

Here is the final script:

```python
import inductiva


machine_group = inductiva.resources.MachineGroup(
    machine_type="c3d-standard-90",
    spot=True,
    data_disk_gb=20)  ## We are adding extra disk space here
machine_group.start()


reef3d = inductiva.simulators.REEF3D()

task = reef3d.run(input_dir="./10_2_3D_Dam_Break_with_Obstacle",
                  on=machine_group,
                  n_vcpus=90)   ## We need to make sure this matches `N 10`

task.wait()
task.download_outputs()

## Make sure we terminate the machine.
machine_group.terminate()
```

Running this script end-to-end should take about 25 minutes, distributed in
the following way:

* about 2 minutes in the preparation stage, including starting the executer VM.
* about 11 minutes of compute time (you can change this by choosing other VMs
or changing the settings of this one.)
* 4 minutes to make the output available for download. This involves zipping
all the data geneeated in executer VM and the move it to your personal storage
space on Inductiva platform, so you can access it at any time.
* about 7 minutes to download the 3.3GB of zipped data. This time depends
on the speed of your network. This is an optional step. You don't need to
download the output data at this point.
* about 3 minutes to unzip the data locally (about 15GB of data when unzipped).
This is another optional step that depends on how fast your local machine is.
* 1 minute to turn off the executor VM. This step is optional but advisable,
since you do not want to be spending money on a machine that is sitting idle.

This is what you are supposed to see:

```
Registering MachineGroup configurations:
> Name:         api-kw1m7e9hs2yxkxy9n2yf4so6r
> Machine Type: c3d-standard-90
> Data disk size:    20 GB
> Number of machines: 1
> Spot:               True
> Estimated cloud cost of machine group: 1.496 $/h
Starting MachineGroup(name="api-kw1m7e9hs2yxkxy9n2yf4so6r"). This may take a few minutes.
Note that stopping this local process will not interrupt the creation of the machine group. Please wait...
Machine Group api-kw1m7e9hs2yxkxy9n2yf4so6r with c3d-standard-90 machines successfully started in 0:00:23.
Task Information:
> ID:                    ggkjuzhivoon56vkozgqxapfk
> Method:                reef3d
> Local input directory: 10_2_3D_Dam_Break_with_Obstacle
> Submitting to the following computational resources:
 >> Machine Group api-kw1m7e9hs2yxkxy9n2yf4so6r with c3d-standard-90 machines
Preparing upload of the local input directory 10_2_3D_Dam_Break_with_Obstacle (399 B).
Input archive size: 640 B
Uploading input archive...
100%|█████████████████████████████████████████████████████████████████████████████████| 640/640 [00:00<00:00, 1.03kB/s]
Local input directory successfully uploaded.
Task ggkjuzhivoon56vkozgqxapfk submitted to the queue of the Machine Group api-kw1m7e9hs2yxkxy9n2yf4so6r with c3d-standard-90 machines.
Simulation metadata logged to: inductiva_output/task_metadata.json
Task ggkjuzhivoon56vkozgqxapfk configurations metadata saved to the tasks metadata file task_metadata.json in the current working directory.
Consider tracking the status of the task via CLI:
	inductiva tasks list --id ggkjuzhivoon56vkozgqxapfk
Or, tracking the logs of the task via CLI:
	inductiva logs ggkjuzhivoon56vkozgqxapfk
Task ggkjuzhivoon56vkozgqxapfk successfully queued and waiting to be picked-up for execution...
Task ggkjuzhivoon56vkozgqxapfk has started and is now running remotely.
Task ggkjuzhivoon56vkozgqxapfk completed successfully.
Downloading simulation outputs to inductiva_output/ggkjuzhivoon56vkozgqxapfk/output.zip.
100%|█████████████████████████████████████████████████████████████████████████████| 3.30G/3.30G [07:11<00:00, 7.64MB/s]
Uncompressing the outputs to inductiva_output/ggkjuzhivoon56vkozgqxapfk.
Terminating MachineGroup(name="api-kw1m7e9hs2yxkxy9n2yf4so6r"). This may take a few minutes.
Machine Group api-kw1m7e9hs2yxkxy9n2yf4so6r with c3d-standard-90 machines successfully terminated in 0:01:11.
```

You can check the stdout of the simulation process in real time by issuing:

```
inductiva logs ggkjuzhivoon56vkozgqxapfk
```

The command line above is also shown in the execution trace, so you can just
copy and paste it to a new terminal (which needs also to have the API key set
as an environment variable). Then, you should see something like this:

```
...
...
7193
simtime: 24.982
timestep: 0.007
Volume 1: 0.235
Volume 2: 1.725
lsmtime: 0.004
.piter: 3  ptime: 0.010
.piter: 3  ptime: 0.010
.piter: 3  ptime: 0.010
umax: 0.238 	 utime: 0.007
vmax: 0.199 	 vtime: 0.007
wmax: 0.277 	 wtime: 0.006
viscmax: 0.000
kinmax: 0.000
epsmax: 0.000
reinitime: 0.005
gctime: 0.002	 average gctime: 0.002
Xtime: 0.013	 average Xtime: 0.014
total time: 626.067449   average time: 0.087
timer per step: 0.073
------------------------------------
7194
simtime: 24.989
timestep: 0.007
Volume 1: 0.235
Volume 2: 1.725
lsmtime: 0.003
.piter: 3  ptime: 0.010
.piter: 3  ptime: 0.010
.piter: 3  ptime: 0.010
umax: 0.237 	 utime: 0.007
vmax: 0.202 	 vtime: 0.006
wmax: 0.276 	 wtime: 0.006
viscmax: 0.000
kinmax: 0.000
epsmax: 0.000
reinitime: 0.005
gctime: 0.002	 average gctime: 0.002
Xtime: 0.011	 average Xtime: 0.014
total time: 626.138361   average time: 0.087
timer per step: 0.071
------------------------------------
7195
simtime: 24.996
timestep: 0.007
Volume 1: 0.235
Volume 2: 1.725
lsmtime: 0.004
.piter: 3  ptime: 0.010
.piter: 3  ptime: 0.009
.piter: 3  ptime: 0.009
umax: 0.235 	 utime: 0.007
vmax: 0.204 	 vtime: 0.006
wmax: 0.276 	 wtime: 0.006
viscmax: 0.000
kinmax: 0.000
epsmax: 0.000
reinitime: 0.005
gctime: 0.002	 average gctime: 0.002
Xtime: 0.013	 average Xtime: 0.014
total time: 626.217546   average time: 0.087
timer per step: 0.079

******************************

modelled time: 25.003

```


## What to read next

If you are interested in Reef3D, you may also be interested in checking the
following related simulators that are also avaiable via Inductiva API:

* [SCHISM](SCHISM.md)
* [SWAN](SWAN.md)
* [SWASH](SWASH.md)
* [XBeach](XBeach.md)