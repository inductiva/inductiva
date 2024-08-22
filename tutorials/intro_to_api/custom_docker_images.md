# Running Simulations with Custom Docker Images

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

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/fds-input-example.zip", unzip=True)

custom_simulator = inductiva.simulators.CustomImage(container_image="docker://inductiva/kutu:fds_v6.8")

task = custom_simulator.run(input_dir=input_dir, commands=["fds mccaffrey.fds"])

task.wait()
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

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/fds-input-example.zip", unzip=True)

custom_simulator = inductiva.simulators.CustomImage(container_image="docker://inductiva/kutu:fds_v6.8")

mpi_config = inductiva.commands.MPIConfig("4.1.6", np=4, use_hwthread_cpus=True)
command = inductiva.commands.Command("fds mccaffrey.fds", mpi_config=mpi_config)

task = custom_simulator.run(input_dir=input_dir, commands=[command])

task.wait()
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

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/fds-input-example.zip", unzip=True)

custom_simulator = inductiva.simulators.CustomImage(container_image="docker://inductiva/kutu:fds_v6.8")

# Running as root is *strongly* discouraged as any mistake (e.g., in
# defining TMPDIR) or bug can result in catastrophic damage to the OS
# file system, leaving your system in an unusable state.
command = inductiva.commands.Command("mpirun --allow-run-as-root -np 4 --use-hwthread-cpus fds mccaffrey.fds")

task = custom_simulator.run(input_dir=input_dir, commands=[command])

task.wait()
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
version, please contact us.