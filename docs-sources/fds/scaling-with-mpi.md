
# Scale FDS with MPI

Fire Dynamics Simulator (FDS) supports parallelism using two methods:
- MPI (Message Passing Interface): Distributes work across multiple meshes by running separate MPI processes.
- OpenMP: Multithreads computations within each mesh.

This tutorial walks you through controlling MPI and OpenMP settings when running FDS simulations on Inductiva, using a small MPI scaling benchmark from the FDS GitHub repository.


## Prerequisites

Before running the simulation, make sure you were able to [run your first FDS simulation](setup-test.md).
You’ll also need to download the input files. You can either:
- **Manually download** them from the [FDS GitHub repository](https://github.com/firemodels/fds/tree/FDS-6.10.1/Validation/MPI_Scaling_Tests/FDS_Input_Files) and place them in a folder named `MpiStrongScalingTest`.
The files starting with `strong_scaling_test` are the ones needed.
**or**
- **Download automatically** [here](https://storage.googleapis.com/inductiva-api-demo-files/fds-tutorials/MpiStrongScalingTest.zip).

## Running a Single-Mesh Simulation

The input folder contains a set of FDS input files representing the same simulation case split into varying number of meshes.
For instance, the `strong_scaling_test_001.fds` has a single mesh, while the `strong_scaling_test_008.fds` is the same case split into 8 meshes, and so forth.

Here’s a script to run the 1-mesh case with one MPI process and one OpenMP thread:

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

> **Note**: The `input_dir` parameter in the `fds.run` method should be set to the path
> to your input files folder.

You can verify the MPI/OpenMP setup by inspecting the stderr logs in the [Inductiva web console](http://console.inductiva.ai) in the task page under Task Logs → stderr, which should show:

```
Number of MPI Processes:  1
Number of OpenMP Threads: 1
```

The task summary looks like this:

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

The simulation took approximately 21 minutes. Let's see how to make it faster using MPI.

## Running with Multiple MPI Processes

To parallelize the simulation using MPI, we'll need:
- A simulation case with multiple meshes (*e.g.*, `strong_scaling_test_008.fds`);
- A machine with at least as many vCPUs as MPI processes. For the 8 mesh case,
we can use the `c2d-standard-8`. You can find more about available machines in the (web console)[https://console.inductiva.ai/machine-groups/instance-types];
- Configure `n_mpi_processes` accordingly.

Here’s an updated script running the 8-mesh case on a machine with 8 vCPUs:

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

The summary for this run:

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

This run took about 6 minutes and 20 seconds — roughly a 4× speedup.
Note that even though we increased parallelism by 8 times, real-world inefficiencies in parallel computing — such as communication and synchronization overhead — limit the actual speed-up.

## Results

Below is a summary of simulation time, cost, and speed-up for the 1-, 8-, and 32-mesh cases. The speed-up is relative to the 1-mesh case.

| MPI Processes | Machine Type     | Time (s) | Cost (USD) | Speed-up   |
|---------------|------------------|----------|------------|------------|
| 1             | c2d-standard-2   | 1275.29  | 0.0057*    | 1.0x       |
| 8             | c2d-standard-8   | 321.59   | 0.0047     | 4.0x       |
| 32            | c2d-standard-32  | 109.35   | 0.0062     | 11.7x      |

By increasing the number of meshes from 8 to 32, we observe a further **2.9× speed-up**.
While the simulation cost does increase slightly, the benefit is you'll **get the results much faster**.

With Inductiva, you have full flexibility to choose the **resources that best match your time and budget constraints.**

> *Note:* The 1-process simulation appears more expensive because it only used 1 of the 2 available vCPUs. You're still billed for the whole machine, even if only one core is active.

## Advanced Configuration

In the previous example, we leveraged hyperthreading by running as many MPI processes as there were virtual CPUs (vCPUs). By default, Google Cloud VMs provide 2 vCPUs per physical CPU core.

```
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-standard-8",
    threads_per_core=1,
    spot=True)
```

Below are the results of running the 1-mesh case with 1 MPI process on a `c2d-standard-2` machine and the 8-mesh case with 8 MPI processes on a `c2d-standard-16` machine:

| MPI Processes | Machine Type     | Time (s) | Cost (USD) | Speed-up   |
|---------------|------------------|----------|------------|------------|
| 1             | c2d-standard-2   | 1295.13  | 0.0058     | 1.0x       |
| 8             | c2d-standard-16  | 199.45   | 0.0057     | 6.5x       |

Compared to the earlier 8-mesh run on a `c2d-standard-8` machine, this configuration achieved a higher speed-up, closer to the theoretical maximum of 8×.
However, due to the use of a more expensive machine, the overall cost of the simulation also increased.

## Key Takeaways

- **FDS uses MPI and OpenMP** to parallelize computations. MPI distributes work across meshes; OpenMP accelerates work within each mesh.
- The number of **MPI processes must not exceed the number of meshes** in the simulation.
- Inductiva gives you **full control** over the number of MPI processes and OpenMP threads — so you can balance performance and cost for your needs.
