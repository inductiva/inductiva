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
    "reef3d-input-example.zip", unzip=True)

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


**control.txt (for DiveMESH):**
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

**ctrl.txt (for Reef3D):**
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

Observe that parameter ```M 10```, which controls the level of parallism, is
set to 4 threads/vCPUs, a very low number. Depending on the number of vCPUs we
effectively wish to use for running the simulation, we will need to manually
change  ```M 10``` on **both files** to match the specs of corresponding the VM. 
We will be using GCP VMs of the c2d family, which tend to provide a very good 
price-performance point, especially when you run them in spot mode. More 
specifically, we will be using the larger ```c2d-highcpu-112``` machines with
112 vCPUs and 2GB of RAM per vCPU (this simulation will not require more than
224GB of RAM). 

We will be setting ```M 10 56``` on both files. Note that this is half the
available number of vCPUs available on the VM (112!). Actually, running on only
half of the vCPUs may provide significantly faster simulation times than running
on all the vCPUs *for some simulators / VM types*. This has to do with the fact
that these VMs run on processors with hyperthreading, running two threads per 
physical core. One of the issues with hyperthreading is that it increases
thread competition for cache (which does is fixed), and this may lead to some
contention issues for I/O heavy simulators, as seems to be the case for Reef3D.

Indeed, Reef3D produces a huge amount of data. As it is currently configured, 
this simulation would produce several dozen gigabytes of data. To reduce that
amount of data produced, we can reduce the rate at which (Paraview) data is
being produced (parameter ```P 30```). Therefore, we will set ```P 30 0.04```
and request Reef3D to "only" generate 25 frames per second (instead of 100
frames per second as it was initially configure). Also, just in case, we will
request our machine to be equipped with a 20GB partition (just for data), 
using the ```data_disk_gb``` parameter of the 
```inductiva.resources.MachineGroup``` class.

Here is the final script:

```python
import inductiva


machine_group = inductiva.resources.MachineGroup(
    machine_type="c2d-highcpu-112",    
    spot=True,
    data_disk_gb=20)
machine_group.start()

reef3d = inductiva.simulators.REEF3D()

task = reef3d.run(
    input_dir="./10_2_3D_Dam_Break_with_Obstacle",
    on=machine_group,
    n_vcpus=56,
    storage_dir="3D_dam_break_with_obstacle")

task.wait()
machine_group.terminate()
task.print_summary()
```

You should see something like this when you run this script end to end (it
should take approximately 20 minutes)

```bash
Registering MachineGroup configurations:
> Name:         api-ktq7n9w81gek0kx82agu3znoj
> Machine Type: c2d-highcpu-112
> Data disk size:    20 GB
> Number of machines: 1
> Spot:               True
> Estimated cloud cost of machine group: 0.814 $/h
 >> You are spending 5.7x less by using spot machines.
Starting MachineGroup(name="api-ktq7n9w81gek0kx82agu3znoj"). This may take a few minutes.
Note that stopping this local process will not interrupt the creation of the machine group. Please wait...
Machine Group api-ktq7n9w81gek0kx82agu3znoj with c2d-highcpu-112 machines successfully started in 0:00:19.
The machine group is using the following quotas:

                      USED BY RESOURCE     CURRENT USAGE     MAX ALLOWED
 cost_per_hour        0.81424              0.81424           20
 total_num_machines   1                    1                 10
 total_num_vcpus      112                  112               150

Task Information:
> ID:                    n90hpbaiaa26jjb1d93zgbqgn
> Method:                reef3d
> Local input directory: 10_2_3D_Dam_Break_with_Obstacle
> Submitting to the following computational resources:
 >> Machine Group api-ktq7n9w81gek0kx82agu3znoj with c2d-highcpu-112 machines
Preparing upload of the local input directory 10_2_3D_Dam_Break_with_Obstacle (392 B).
Input archive size: 637 B
Uploading input archive...
100%|███████████████████████████████████████████████████████████████████████████████| 637/637 [00:00<00:00, 1.37kB/s]
Local input directory successfully uploaded.
Task n90hpbaiaa26jjb1d93zgbqgn submitted to the queue of the Machine Group api-ktq7n9w81gek0kx82agu3znoj with c2d-highcpu-112 machines.
Number of tasks ahead in the queue: 0
Simulation metadata logged to: inductiva_output/task_metadata.json
Task n90hpbaiaa26jjb1d93zgbqgn configurations metadata saved to the tasks metadata file task_metadata.json in the current working directory.
Consider tracking the status of the task via CLI:
	inductiva tasks list --task-id n90hpbaiaa26jjb1d93zgbqgn
Or, tracking the logs of the task via CLI:
	inductiva logs n90hpbaiaa26jjb1d93zgbqgn
Task n90hpbaiaa26jjb1d93zgbqgn successfully queued and waiting to be picked-up for execution...
The task n90hpbaiaa26jjb1d93zgbqgn is about to start.
Task n90hpbaiaa26jjb1d93zgbqgn has started and is now running remotely.
 Task n90hpbaiaa26jjb1d93zgbqgn completed successfully.
Downloading stdout and stderr files to inductiva_output/n90hpbaiaa26jjb1d93zgbqgn...
Partial download completed to inductiva_output/n90hpbaiaa26jjb1d93zgbqgn.
Successfully requested termination of MachineGroup(name="api-ktq7n9w81gek0kx82agu3znoj").
Termination of the machine group freed the following quotas:

                      FREED BY RESOURCE     CURRENT USAGE     MAX ALLOWED
 cost_per_hour        0.81424               0                 20
 total_num_machines   1                     0                 10
 total_num_vcpus      112                   0                 150


Task status: success
Wall clock time:  0:13:51
Time breakdown:
	Input upload:              0.63 s
	Time in queue:             24.34 s
	Container image download:  3.37 s
	Input download:            0.08 s
	Input decompression:       0.00 s
	Computation:               0:07:38
	Output compression:        0:05:04
	Output upload:             40.08 s
Data:
	Size of zipped output:    3.09 GB
	Size of unzipped output:  7.51 GB
	Number of output files:   35751
```

While the script is running, you can check the stdout of the simulation process
in real time by issuing (please change to the id of your task):

```bash
inductiva logs n90hpbaiaa26jjb1d93zgbqgn
```

The command line above is shown in the execution trace, so you can just
copy and paste it to a new terminal (which needs also to have the API key set
as an environment variable).

Once the script finishes, and you see the summary of the task, you can download
the resulting files. Observe in the summary above that quite some data was 
produced: 35751 files for a total of 7.51 GB after uncompressed.

You can download the (zipped) data by creating a simple script such as this:
```python
import inductiva

# You can retreive a Task by ID
task = inductiva.tasks.Task("n90hpbaiaa26jjb1d93zgbqgn")

task.download_outputs()
```

Or, perhaps more conveniently, you can use Inductiva's Command Line Interface
(CLI) to list the contents of your personal remote storage and then download
them:
```bash
inductiva storage list
```
which should get you something like this (you may have other contents)

```
 NAME                          SIZE      CREATION TIME
 3D_dam_break_with_obstacle/   3.09 GB   26 Jun, 07:39:11
```

And you can donwload the folder produced by this task by doing:

```bash
inductiva tasks download n90hpbaiaa26jjb1d93zgbqgn
```
## What to read next

If you are interested in Reef3D, you may also be interested in checking the
following related simulators that are also avaiable via Inductiva API:

* [SCHISM](SCHISM.md)
* [SWAN](SWAN.md)
* [SWASH](SWASH.md)
* [XBeach](XBeach.md)