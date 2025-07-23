# Scale your FDS simulations with MPI

There are two ways of exploring multi-processing and parallelism in FDS simulations: MPI and OpenMP. The number of MPI processes
can be used to parallelize computation for several meshes (*e.g.*, if your problem is decomposed in 4 meshes, you can leverage 4 MPI processes to speed up your simulation).

As such, the number of MPI processes that can be used in a given simulation is limited by the number of meshes in the problem.
OpenMP is used to parallelize computations within a given mesh.

While running your FDS via Inductiva, you retain full control of the number of MPI processes and OpenMP threads used by the simulation by settings two parameters when submitting your simulation.
In this tutorial, we'll walk through how to configure those parameters using a small MPI scaling benchmark obtained from the [FDS GitHub repository](https://github.com/firemodels/fds/tree/FDS-6.10.1/Validation/MPI_Scaling_Tests).


## Prerequisites

Before running the simulation, you should have been able to [run your first FDS simulation](setup-test.md).
You’ll also need to download the required input files. You can either:

- **Manually download** them from the [FDS GitHub repository](https://github.com/firemodels/fds/tree/FDS-6.10.1/Validation/MPI_Scaling_Tests/FDS_Input_Files) and place them in a folder named `MpiStrongScalingTest`.
The files starting with `strong_scaling_test` are the ones needed.
**or**
- **Download automatically** using the link provided [here](https://storage.googleapis.com/inductiva-api-demo-files/fds-tutorials/MpiStrongScalingTest.zip).

## Running the single mesh case

The input folder includes a series of FDS inputs files, with a simple simulation divided in a different number of meshes.
For instance, the `strong_scaling_test_001.fds` has a single mesh, while the `strong_scaling_test_008.fds` is the same case divided in 8 meshes.

The following script can be used to run the simulation:

```python
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-standard-2",
    spot=True)

# Initialize the Simulator
fds = inductiva.simulators.FDS( \
    version="6.9.1")

# Run simulation
task = fds.run( \
    input_dir="/Path/to/MpiStrongScalingTest",
    sim_config_filename="strong_scaling_test_001.fds",
    n_mpi_processes=1,
    n_omp_threads=1,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

> **Note**: The `input_dir` parameter in the `fds.run` method should be set with the correct path
> to the folder containing the input files.

Since the `strong_scaling_test_001.fds` refers to single mesh case, we use the `n_mpi_processes` parameter to configure the simulation to run with a single MPI process.
As such, with the above script we run a single mesh case with one MPI process a GCP `c2d-standard-2` machine (a machine with 2 vCPUs). Note that each MPI process is using
a single OpenMP thread.

This can be confirmed by inspecting the stderr logs of the FDS command. Navigate to the [web console](http://console.inductiva.ai) and find the task you just ran.
Navigate to the `Task Logs` tab and select stderr. You should see the following lines:

```
Number of MPI Processes:  1
Number of OpenMP Threads: 1
```

The outputs of `task.print_summary()` are the following:

```
Task status: Success

Timeline:
        Waiting for Input         at 23/07, 16:55:18      4.989 s
        In Queue                  at 23/07, 16:55:23      42.639 s
        Preparing to Compute      at 23/07, 16:56:06      5.442 s
        In Progress               at 23/07, 16:56:11      1275.291 s
                └> 1275.151 s      /opt/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux strong_scaling_test_001.fds
        Finalizing                at 23/07, 17:17:26      0.538 s
        Success                   at 23/07, 17:17:27

Data:
        Size of zipped output:    45.84 KB
        Size of unzipped output:  1.48 MB
        Number of output files:   19

Estimated computation cost (US$): 0.0057 US$
```

The simulation took around 21 minutes to run. Let's see how to make it faster using MPI.

## Running with multiple MPI processes

To parallelize the simulation using MPI, we'll need:
- A simulation case with multiple meshes;
- A machine with multiple vCPUs;
- Configure the `n_mpi_processes` parameter

To run the 8 mesh case (`strong_scaling_test_008.fds`), we need a machine with at least 8 vCPUs,
such as the `c2d-standard-8`. You can find more about available machines in the (web console)[https://console.inductiva.ai/machine-groups/instance-types].

The updated script is as follows:

```python
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-standard-8",
    spot=True)

# Initialize the Simulator
fds = inductiva.simulators.FDS( \
    version="6.9.1")

# Run simulation
task = fds.run( \
    input_dir="/Path/to/MpiStrongScalingTest",
    sim_config_filename="strong_scaling_test_008.fds",
    n_mpi_processes=8,
    n_omp_threads=1,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

The task summary:

```
Task status: Success

Timeline:
        Waiting for Input         at 23/07, 17:03:39      0.803 s
        In Queue                  at 23/07, 17:03:39      51.276 s
        Preparing to Compute      at 23/07, 17:04:31      4.398 s
        In Progress               at 23/07, 17:04:35      321.586 s
                └> 321.427 s       /opt/openmpi/4.1.6/bin/mpirun --np 8 --use-hwthread-cpus /opt/fds/Build/ompi_gnu_linux/fds_ompi_gnu_linux strong_scaling_test_008.fds
        Finalizing                at 23/07, 17:09:57      0.429 s
        Success                   at 23/07, 17:09:57

Data:
        Size of zipped output:    66.40 KB
        Size of unzipped output:  1.57 MB
        Number of output files:   75

Estimated computation cost (US$): 0.0047 US$
```

The simulation now ran in about 6 minutes and 20s. Around 4x faster. Note that, even though we increase parallelism by 8 times, in practice,
inefficiences related with communication and sincronization don't allow for such speed ups.

## Results

Here we present the simulation time and cost of the above examples, as well as the 32 mesh case.

| MPI proc | Machine Type    | Time (s) | Cost (US$) |
|----------|-----------------|----------|------------|
| 1        | c2d-standard-2  | 1275.291 | 0.0057     |
| 8        | c2d-standard-8  | 321.586  | 0.0047     |
| 32       | c2d-standard-32 | 109.349  | 0.0062     |

