# Run GRIR443 on an MPI
QuantumESPRESSO is designed from the ground up with parallelism in mind, enabling it to efficiently distribute both data and 
computational workload across multiple processors and nodes. This is particularly beneficial for demanding tasks such as plane-wave 
DFT (Density Functional Theory) calculations, which involve large Hamiltonian matrices and extensive Fourier transforms. By 
leveraging **MPI (Message Passing Interface)**, QuantumESPRESSO can parallelize key components like k-point sampling, band 
structure calculations, and FFTs, significantly reducing runtime and enabling simulations of larger and more complex systems.

With Inductiva, you can easily distribute the computational load across multiple nodes by making minor modifications to the single-node scripts we presented earlier in this [tutorial](run-grir443-benchmark).

All you need to do is update the resource allocation from `MachineGroup` to `MPICluster`, as follows:

```diff
# Allocate a multi-machine MPI cluster on Google Cloud Platform
- cloud_machine = inductiva.resources.MachineGroup(
+ cloud_machine = inductiva.resources.MPICluster(
    machine_type="c3d-highcpu-360",
+   num_machines=2,
    spot=True)
```

In this tutorial, we’ll demonstrate how leveraging an MPI cluster — both with and without hyperthreading enabled — can impact 
the performance of the [GRIR443 benchmark](https://github.com/QEF/benchmarks/tree/master/GRIR443).

## Running the Simulation
The script for running the GRIR443 simulation on a multi-node MPI cluster remains largely the same as the single-machine version. The c
ore structure, including defining commands, initializing the simulator, and launching the job, does not change. The main differences lie 
in resource allocation and MPI configuration, where we adjust the number of machines and the number of MPI processes (threads) to match 
the parallel environment.

In this example, we request two `c3d-highcpu-360` machines, each with 360 vCPUs, connected as an MPI cluster. This setup provides a total 
of 720 threads (2 × 360).

To reflect this, we update the `MPIConfig` object by setting the np parameter to 720:

```python 
mpi_config = MPIConfig(
    version="4.1.6",
    np=720,
    use_hwthread_cpus=True
)
```

Note that `use_hwthread_cpus` is set to `True`. While this is the default and does not strictly need to be specified, we include it here for clarity, as its purpose will become clearer later in the tutorial.

Building on the single-machine setup, here's how to launch the GRIR443 simulation on a 2-machine MPI cluster
using the Inductiva Python API:

```python
import inductiva
from inductiva.commands import MPIConfig, Command

# Allocate a multi-machine MPI cluster on Google Cloud Platform
cloud_machine = inductiva.resources.MPICluster(
    machine_type="c3d-highcpu-360",
    num_machines=2,
    spot=True)

mpi_config = MPIConfig( \
    version="4.1.6",
    np=720,
    use_hwthread_cpus=True)

# List of commands to run
commands = [
    Command("pw.x -i grir443.in", mpi_config=mpi_config),
]

# Initialize the Simulator
qe = inductiva.simulators.QuantumEspresso(\
    version="7.4.1")

# Run simulation
task = qe.run( \
    input_dir="GRIR443",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()
task.download_outputs()
task.print_summary()
```

This simulation runs in approximately **415** seconds using a two-machine MPI cluster.

By comparison, running it on a single machine took around 800 seconds, resulting in a **1.93× speedup**.

## Scaling to Larger MPI Clusters
A major benefit of using Inductiva is how easily you can scale your simulations to larger MPI clusters with minimal code changes. This is achieved simply by updating the `num_machines` parameter in the MPICluster configuration and the `np` value in the MPIConfig object.

Below are the results of running this simulation on increasingly larger MPI clusters, compared to the single-node setup (shown in the first row):

| Machines (`num_machines`) | MPI Processes (`np`) | Time (s) | Speedup | Cost (USD) |
|---------------------------|----------------------|----------|---------|------------|
| 1                         | 360                  | 800      | —       | 0.76       |
| 2                         | 720                  | 415      | 1.93×   | 0.81       |
| 4                         | 1440                 | 275      | 2.91×   | 1.11       |
| 8                         | 2880                 | 188      | 4.26×   | 1.50       |

> **Note**: Hyperthreading is enabled.

We achieve significant speed-ups by changing just a couple of lines of code. All this improved performance comes without a dramatic increase in the simulation cost. It’s impressive that you can run your simulations more than **four times faster** by roughly doubling the cost. If saving time is a priority, this is an excellent trade-off!

However, as expected, we observe diminishing returns when increasing the number of nodes. One pattern we have consistently seen — including with other simulators — is that scaling Inductiva’s MPICluster up to 8 machines generally improves performance, but beyond that point, it drops off sharply.

There are two main reasons for this. The primary and most fundamental reason relates to the internode connections within the cloud infrastructure supporting Inductiva. Unlike traditional HPC systems, the cloud was not designed for such a high level of connectivity. Moreover, since cloud resources are shared among many users, internode connections can become unpredictably busy. As a result, the more nodes in an MPI cluster, the greater the chance that one slow connection will bottleneck the entire computation. Unfortunately, there is little we can do about this except hope for favorable internode performance.

The second reason concerns the default settings of the Google Cloud VMs used in these tests. By default, Google VMs have hyperthreading enabled, meaning two vCPUs (threads) run per physical core. At very high thread counts per node, this can lead to memory access bottlenecks, further slowing down the simulation. This is one reason traditional HPC infrastructures often disable hyperthreading. Fortunately, Inductiva allows you to disable hyperthreading if desired. Let’s explore how to do that and rerun the tests.

## Turning Off Hyperthreading
Repeating the tests above without hyperthreading requires only two small changes. First, when initializing the `MPICluster`, we need to explicitly specify that we want only one thread per core. This is done by setting the `threads_per_core` parameter:

```python
cloud_machine = inductiva.resources.MPICluster(
    machine_type="c3d-highcpu-360",
    num_machines=2,
    threads_per_core=1,
    spot=True)
```

Next, we update the `MPIConfig` to reflect this change. There are two important adjustments:
- Set `np` to half the previous value because the number of available threads is now halved.
- Set `use_hwthread_cpus` to `False` since each physical core will now support only a single thread/vCPU.

For example, with two `c3d-highcpu-360` nodes (which have 2 × 180 = 360 physical cores), the configuration becomes:

```python
mpi_config = MPIConfig(
    version="4.1.6",
    np=360,
    use_hwthread_cpus=False)
```

**Does this improve performance?** Let’s see the results:

| Machines (`num_machines`) | MPI Processes (`np`) | Cores | Time (s) | Speedup | Cost (USD) |
|---------------------------|----------------------|-------|----------|---------|------------|
| 1                         | 180                  | 180   | 787      | —       | 0.75       |
| 2                         | 360                  | 360   | 399      | 1.97×   | 0.77       |
| 4                         | 720                  | 720   | 210      | 3.74×   | 0.83       |
| 8                         | 1440                 | 1440  | 154      | 5.11×   | 1.24       |

The results are clear: all runs *without* hyperthreading are **faster** than their counterparts with hyperthreading enabled — in some 
cases by **more than 20%**. The speedups are also consistently better, achieved at **very reasonable costs**. In fact, you can 
achieve speedups of over 5× by increasing the cost by only about 65%!

**How does this performance compare to a proper supercomputer?** If you’re curious, find out more [here](benchmarks) to see how these numbers compare with those from Fugaku, one of the fastest supercomputers on Earth.

