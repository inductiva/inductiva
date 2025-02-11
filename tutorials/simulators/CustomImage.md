# Custom Docker Images

In this tutorial, we will guide you through running simulations with custom
Docker images. This can be particularly useful if you need to use a specific
version of software that isn't available on our Docker Hub. Using this feature,
you can run your own Docker image on the platform.

This tutorial will focus on two main points:
- The simulator you want to use
- The inputs required by that simulator

## CustomImage Simulator

The `CustomImage` simulator allows you to specify a `container_image` pointing
to your custom Docker image. Besides the `container_image`, this simulator
accepts a list of commands to run inside the container. These are the only
special parameters for this simulator. The rest are the same as other simulators.

Here’s an example of how to use the `CustomImage` simulator:

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c3d-standard-90")
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/fds-input-example.zip", unzip=True)

custom_simulator = inductiva.simulators.CustomImage(container_image="docker://inductiva/kutu:fds_v6.8")

task = custom_simulator.run(input_dir=input_dir, commands=["fds mccaffrey.fds"],
                            on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()
```

This basic example demonstrates how to run a custom Docker image. To leverage
this feature fully, we introduce the `Command` and `MPIConfig` classes, which
enable you to run commands on MPI clusters supporting multiple MPI versions.

## Command and MPIConfig

In the previous example, `commands` is a list of strings to be run inside the
container. For more advanced usage, such as running a command with a specific
MPI version, use the `Command` and `MPIConfig` classes. A list of `commands` can
include `Command` instances, each potentially having an `MPIConfig`. Here’s an
example:

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/fds-input-example.zip", unzip=True)

custom_simulator = inductiva.simulators.CustomImage(container_image="docker://inductiva/kutu:fds_v6.8")

mpi_config = inductiva.commands.MPIConfig("4.1.6", np=4, use_hwthread_cpus=True)
command = inductiva.commands.Command("fds mccaffrey.fds", mpi_config=mpi_config)

task = custom_simulator.run(input_dir=input_dir, commands=[command],
                            on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()
```

This example runs four instances of your image with the command `fds mccaffrey.fds`
using MPI version 4.1.6, equivalent to:

```sh
mpirun -np 4 -use_hwthread_cpus apptainer run ... fds mccaffrey.fds
```

Alternatively, you can run MPI directly inside your container with:

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/fds-input-example.zip", unzip=True)

custom_simulator = inductiva.simulators.CustomImage(container_image="docker://inductiva/kutu:fds_v6.8")

command = inductiva.commands.Command("mpirun -np 4 --use-hwthread-cpus fds mccaffrey.fds")

task = custom_simulator.run(input_dir=input_dir, commands=[command],
                            on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()
```

This runs:

```sh
apptainer run ... mpirun -np 4 -use_hwthread_cpus fds mccaffrey.fds
```

Ensure the MPI version in your container matches the one specified in `MPIConfig`.
For instance, if your container has OpenMPI 1.10.6, choose a compatible OpenMPI
version (like 1.10.7) in `MPIConfig`.

Currently, we support OpenMPI versions 1.10.7 and 4.1.6. If you need a different
version, please contact us at [contact@inductiva.ai](mailto:contact@inductiva.ai).

## Advanced Tutorial: Running our Reef3D Advanced example

Let's now run a more advanced example. For this we will be using our 
[advanced tutorial on Reef3D](Reef3D.md#advanced-tutorial-3d-dam-break-scenario-with-obstacle)
but we will be running it using our `CustomImage` simulator.

We will be skipping some explenations related with Reef3D and will focus more on
the parts related with the `CustomImage` simulator.


### Aim of This Tutorial

This tutorial will guide you through running a simulation on our platform using
a custom Docker image so you can run your own docker images in the future.

### Prerequisites

1. **Download Input Files**: Get the input files from the
[Reef3D tutorials](https://github.com/REEF3D/REEF3D/tree/master/Tutorials/REEF3D_CFD/10_2%203D%20Dam%20Break%20with%20Obstacle).

	**Directory Structure**:
   ```bash
   ls -lasgo .
    total 16
    0 drwxrwxr-x@  4   128 Sep  4 08:46 .
    0 drwxrwxr-x@ 19   608 Nov  5 09:03 ..
    8 -rw-rw-r--@  1   142 Sep  4 08:46 control.txt
    8 -rw-rw-r--@  1   141 Sep  4 08:46 ctrl.txt
   ```

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

	M 10 4    ---- defines the nr. of processors for parallel computations (4)
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
	N 41 25.0    ---- set the maximum modeled time (25 seconds).
	N 45 50000
	N 47 0.2
	M 10 4    ---- defines the nr. of processors for parallel computations (4)
	P 10 1
	P 30 0.01    ---- defines the rate of paraview results (1 frame per 0.01 s)
	T 10 0
	W 22 -9.81
	```

### Overview

Here’s the code you'll be working on as we progress through the tutorial. Don’t
worry if it doesn’t all make sense right now; everything will become clearer
in the upcoming steps.

```python
import inductiva

machine_group = inductiva.resources.MachineGroup(
					machine_type="c2d-highcpu-112",
					spot=True,
					data_disk_gb=20,
					auto_resize_disk_max_gb=250)

machine_group.start()

input_dir = "path/to/10_2 3D Dam Break with Obstacle"

#Choose your simulator
customImage = inductiva.simulators.CustomImage(
    container_image="docker://inductiva/kutu:reef3d_v24.02")

# Define the MPI configuration
mpi_config = inductiva.commands.MPIConfig("4.1.6", np=56, use_hwthread_cpus=False)

# Define the commands to run with the MPI configuration
command = inductiva.commands.Command("/REEF3D/bin/REEF3D", mpi_config=mpi_config)

# List of commands to run our simulation
commands = [
    "/DIVEMesh/bin/DiveMESH",
    command
]

task = customImage.run(
	input_dir=input_dir,
    commands=commands,
	on=machine_group,
	n_vcpus=56,
	use_hwthread=False,
	storage_dir="3D_dam_break_with_obstacle")

task.wait()
machine_group.terminate()

task.print_summary()
```

### Step 1: Adjust Simulation Parameters
For a faster simulation, modify the following parameters in both files:

- **Level of parallelism (`M 10`)**: 56
- **Simulation time (`N 41`)**: 25.0
- **Paraview results rate (`P 30`)**: 0.01


### Step 2: Running the Simulation

#### a. Configure and Start Machine

1. **Pick your machine**:

    Our `CustomImage` simulator will work with any machine type. So, we decided
    to use the same machine as in the previous tutorial, the `c2d-highcpu-112`
    with the same configuration.

	```python
	import inductiva
	machine_group = inductiva.resources.MachineGroup(
						machine_type="c2d-highcpu-112",
						spot=True,
						data_disk_gb=20,
						auto_resize_disk_max_gb=250)
	```
	**Note**: `spot` machines are a lot cheaper but can be terminated by the
	provider if needed.

2. **Start your machine**
	```python
	machine_group.start()
	```

#### b. Define Simulation Inputs

1. **Specify Simulation Directory**:
	Let's start by defining a variable that points to the `10_2_3D_Dam_Break_with_Obstacle`
	folder where all your simulation files are located.

	```python
	input_dir = "./10_2_3D_Dam_Break_with_Obstacle"
	```

2. **Choose Your Simulator**:

This step is straightforward. We will be using the `CustomImage` simulator to
run our custom Docker image. The `container_image` parameter points to the
Docker image we want to use.

```python
customImage = inductiva.simulators.CustomImage(
    container_image="docker://inductiva/kutu:reef3d_v24.02")
```

3. **Define Commands**:
    Our `CustomImage` simulator accepts a list of commands to run inside the
    container. We will be running two commands: `DiveMESH` and `REEF3D`, present
    in the Reef3D image.

    ```python
    # Define the MPI configuration
    mpi_config = inductiva.commands.MPIConfig("4.1.6", np=56, use_hwthread_cpus=False)

    # Define the commands to run with the MPI configuration
    command = inductiva.commands.Command("/REEF3D/bin/REEF3D", mpi_config=mpi_config)

    # List of commands to run our simulation
    commands = [
        "/DIVEMesh/bin/DiveMESH",
        command
    ]
    ```

    Lets talk about the `MPIConfig` for a command. The `MPIConfig` class allows you
    to run your command with a specific MPI version. In this case, we are using
    OpenMPI version 4.1.6 with 56 processes. Creating your command with an `MPIConfig`
    will parallelize said command across the specified number of processes.

    **Note:** Not all commands can be parallelized. Make sure the command you are
    running supports parallelization.

#### c. Run Your Simulation

1. **Run the simulation**:
	We now have everything we need to run our simulation.
	```python
	task = customImage.run(
        input_dir=input_dir,
        commands=commands,
        on=machine_group,
        n_vcpus=56,
        use_hwthread=False,
        storage_dir="3D_dam_break_with_obstacle")
	```

	As you can see, this is very similar to the Reef3D example. The only
	difference here is that we need to specify the specific commands required
	to run the simulation.

2. **Wait for the simulation to finish**:
	That is it. Our simulation is now running on the cloud. Every step from now
    on is the same as the Reef3D example.

	```python
	task.wait()
	```
	**Note**: run `inductiva logs task_id` to check the `stdout` of the simulation
	process in real time.

3. **Terminate Machine**:
	Once our simulation completes, we can/should terminate our machine to save on costs.
	If you forget, don’t worry—we’ve got you covered. By default, the machine will 
	automatically shut down if idle for 30 minutes with no simulation running.


	```python
	machine_group.terminate()
	```

4. **Check your simulation summary**:
	Now that our simulation is complete, we can print a summary with key details, 
	including execution times, generated outputs, and more.

	```python
	task.print_summary()
	```

### Conclusion

Congratulations! You've successfully completed the advanced Reef3D example using
a custom Docker image on our platform. Through this tutorial, you’ve gained
practical experience in configuring and running simulations with custom container
images, even setting up parallel computations with MPI. 

With these skills, you're now prepared to run your own custom Docker images in
future simulations, tailoring the setup to suit different types of computational
requirements. Don’t hesitate to explore other simulations and modify parameters to deepen your understanding of both Reef3D and our custom Docker-based simulator setup. 

Happy simulations!
