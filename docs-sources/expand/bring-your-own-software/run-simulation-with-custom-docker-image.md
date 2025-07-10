# Run Simulation with Custom Docker Image
This guide demonstrates how to run a simulation using a custom Docker image.

By the end, you'll be able to use the `CustomImage` simulator using your own public or private image.

## The `CustomImage` Simulator
The `CustomImage` enables you to run simulations using any Docker image of your choice by specifying the `container_image` parameter. 
This gives you full control over the simulation environment.

In addition to `container_image`, you can provide a list of commands to execute within the container. These are the only parameters unique 
to the `CustomImage` simulator; all other configuration options work just like those for the simulators integrated into our platform.

To illustrate, we’ll use a Docker image of the Fire Dynamics Simulator (FDS), which is publicly available in our Docker Hub repository,
[Kutu](https://hub.docker.com/r/inductiva/kutu). Below is an example of how to run a simulation using the `CustomImage` simulator:

```python
import inductiva

# Allocate cloud machine
machine_group = inductiva.resources.MachineGroup("c2d-highcpu-4")
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/fds-input-example.zip", unzip=True)

custom_simulator = inductiva.simulators.CustomImage(container_image="docker://inductiva/kutu:fds_v6.8")

task = custom_simulator.run(input_dir=input_dir, commands=["fds mccaffrey.fds"],
                            on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
```

To adapt this script for your own custom image, simply update the `container_image` and adjust the `commands` as needed for your simulation.

## Commands and MPI

### Running a Command with a specific MPI version
A list of `commands` can include `Command` instances, each optionally associated with an `MPIConfig`. Below is an example demonstrating this:

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
task.download_outputs()

machine_group.terminate()
```

This example runs four instances of your image using the command `fds mccaffrey.fds`, with MPI version 4.1.6. This is equivalent to:

```
mpirun -np 4 -use_hwthread_cpus apptainer run ... fds mccaffrey.fds
```


### Running MPI directly inside the Container
Alternatively, you can run MPI directly within your container as shown below:

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
task.download_outputs()

machine_group.terminate()
```

This is equivalent to:

```sh
apptainer run ... mpirun -np 4 -use_hwthread_cpus fds mccaffrey.fds
```


### Important Notes

- **MPI Version Compatibility**: Ensure that the MPI version in your container aligns with the version specified in `MPIConfig`. For example, if your container
  has OpenMPI 1.10.6, select a compatible version (such as 1.10.7) in `MPIConfig`.
- **Supported MPI Versions**: We currently support OpenMPI versions 1.10.7 and 4.1.6. If you require a different version, please [Contact Us](mailto:support@inductiva.ai).
