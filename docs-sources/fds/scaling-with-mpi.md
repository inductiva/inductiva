# Scale FDS with MPI
Fire Dynamics Simulator (FDS) supports parallelism using two methods:
- **MPI (Message Passing Interface)**: Distributes the workload across multiple meshes by running separate MPI processes.
- **OpenMP**: Enables multi-threading within each mesh.

In this tutorial, you'll learn how to configure configure MPI settings for running FDS simulations on Inductiva, using a benchmark case from the official FDS GitHub repository designed to test MPI scaling.

## Prerequisites
Before proceeding, make sure you've successfully run your [first FDS simulation](quick-start.md).

Next, download the required input files. You can do this in one of two ways:
- **Manually download** the files from the [FDS GitHub repository](https://github.com/firemodels/fds/tree/FDS-6.10.1/Validation/MPI_Scaling_Tests/FDS_Input_Files) and place them in a folder named `MpiStrongScalingTest`. Only the files that start with `strong_scaling_test` are needed.
**or**
- **Download automatically** [here](https://storage.googleapis.com/inductiva-api-demo-files/fds-tutorials/MpiStrongScalingTest.zip).

## Running a Single-Mesh Simulation
The `MpiStrongScalingTest` folder contains a set of FDS input files representing the same simulation case split into varying number of meshes.
For instance, the `strong_scaling_test_001.fds` has a single mesh, while the `strong_scaling_test_008.fds` is the same case split into 8 meshes, and so forth.

Let’s start by running the **1-mesh** case using **1 MPI process** and **1 OpenMP thread**.

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

> **Note**: Make sure the `input_dir` parameter is set to the path of your `MpiStrongScalingTest` folder.

To verify the MPI/OpenMP setup, check the stderr logs in the [Inductiva Console](http://console.inductiva.ai) under the task details (**Task Logs → stderr.txt**), which should show:

```
Number of MPI Processes:  1
Number of OpenMP Threads: 1
```

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown below.

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

Estimated Task Compute Cost = 0.0057 US$
Task Orchestration Fee = 0.01 US$
Total Estimated Cost = 0.0157 US$
Learn more about costs at: https://inductiva.ai/guides/how-it-works/basics/how-much-does-it-cost
```

The simulation took approximately 21 minutes. Let's use MPI to make it faster.

## Running with Multiple MPI Processes
To run in parallel using MPI:
- Use an input file with multiple meshes (e.g., `strong_scaling_test_008.fds`)
- Select a machine with at least as many virtual CPUs (vCPUs) as the number of MPI processes (for the 8-mesh case, `c2d-standard-8` is a good fit)
- Set `n_mpi_processes` parameter to match the number of meshes

Below is an example running the 8-mesh case on an 8-vCPU machine:

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

> ⚠️ **Note on vCPUs and Hyperthreading**: In most cloud environments (like Google Cloud), a vCPU corresponds to a thread, not a full physical core. So a `c2d-standard-8` machine has 8 vCPUs, which typically means 4 physical cores with hyperthreading enabled.

## Results
The table below compares performance across different MPI configurations for the 1-, 8-, and 32-mesh cases. Speed-up is calculated relative to the 1-mesh baseline.

| Machine Type     | MPI Processes | Execution Time  | Estimated Cost (USD) | Speed-up |
|------------------|---------------|-----------------|----------------------|----------|
| c2d-standard-2   | 1             | 21 minutes, 15s | 0.0057               | 1.0×     |
| c2d-standard-8   | 8             | 5 minutes, 22s  | 0.0047               | 4.0×     |
| c2d-standard-32  | 32            | 1 minute, 49s   | 0.0062               | 11.7×    |

Increasing parallelism from 1 to 8 MPI processes results in a **4× speed-up**, reducing runtime from over 21 minutes to just over 5 minutes. This demonstrates a clear performance gain. However, due to typical overheads in parallel computing, such as communication and synchronization, the speed-up is not perfectly linear.

Scaling further from 8 to 32 processes yields an additional 2.9× speed-up, cutting total runtime to under 2 minutes — **over 11× faster** than the baseline.

While the simulation cost increases slightly with more resources, the main benefit is that you’ll get results much faster.

With Inductiva, you have full flexibility to choose the computational resources that best match your time and budget constraints.

> **Note**: The 1-process simulation appears more expensive because it used only 1 of the 2 available vCPUs. You’re billed for the full machine, regardless of how many cores are utilized.

## Advanced Setup: Disabling Hyper-threading
In the previous examples, we assigned one MPI process per virtual CPU (vCPU), which means the simulations ran on hyperthreads, not physical CPU cores.

However, in traditional HPC environments, it's common practice to run one MPI process per physical core, with hyper-threading disabled. This avoids resource contention and can lead to more predictable and consistent performance.

By default, Google Cloud VMs provide 2 vCPUs per physical core, so hyper-threading is enabled. To disable hyper-threading and use only physical cores, configure the machine group with `threads_per_core=1`:

```
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-standard-16",
    threads_per_core=1,
    spot=True)
```

Here are the results of running without hyper-threading the 1-mesh case with 1 MPI process on a `c2d-standard-2` machine and the 8-mesh case with 8 MPI processes on a `c2d-standard-16` machine:

| Machine Type     | MPI Processes | Execution Time | Estimated Cost (USD) | Speed-up |
|------------------|---------------|----------------|----------------------|----------|
| c2d-standard-2   | 1             | 21 min, 35 s   | 0.0058               | 1.0x     |
| c2d-standard-16  | 8             | 3 min, 19 s    | 0.0057               | 6.5x     |

Compared to the earlier 8-mesh run on a `c2d-standard-8` machine (with hyper-threading enabled), this configuration achieved a higher speed-up, closer to the theoretical maximum of 8×. Despite using a more expensive machine type, the overall cost remained similar, making this setup more efficient in terms of time-to-solution.

## Key Takeaways

- **FDS supports parallelization through MPI and OpenMP**: MPI distributes the workload across multiple meshes, while OpenMP accelerates computations within each mesh.
- The number of MPI processes **must not exceed the number of meshes** defined in your simulation.
- **Running one MPI process per vCPU** (i.e., with hyper-threading enabled) is simple to configure and cost-effective for light to moderate workloads.
- **Disabling hyper-threading** allows MPI processes to run on physical cores, which is common in traditional HPC environments. This setup reduces resource contention and often improves performance consistency.
- Inductiva gives you **full control over MPI and OpenMP settings**, allowing you to optimize for performance, cost, or a balance of both.