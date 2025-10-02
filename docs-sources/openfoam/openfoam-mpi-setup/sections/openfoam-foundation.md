# OpenFOAM-Foundation
OpenFOAM-Foundation runs MPI without extra flags, meaning MPI enforces a strict rule: **the number of processes requested must be less than or equal to the number of physical CPU cores available**.

For example, on a `c2d-highcpu-4` machine, which has 2 physical cores and 4 vCPUs, MPI allows at most 2 processes. This restriction means your simulation domain can only be decomposed into a number of parts equal to or fewer than the physical cores on your computational resource.

## Understanding the Default Behaviour
Let’s revisit the tutorial case from [here](../quick-start.md), running the exact same simulation but using default computational resource settings:

```python
import inductiva

cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True)
```

The simulation is originally split into **6 sub-domains**, so CPU utilization peaks around **37%**:

![CPU Usage](../_static/foundation_6_vcpus.png)

This occurs because the computational resources are configured with `threads_per_core=2` (hyper-threading enabled), which is the default setting for virtual machines (learn more [here](https://inductiva.ai/guides/how-it-works/machines/hyperthreading)). Although the machine has 16 vCPUs, the simulation only uses 6 vCPUs, which explains the low CPU utilization.

You can maximize CPU usage in two ways:
1. Use all available vCPUs (default mode: hyper-threading enabled, with 2 vCPUs per physical core)
2. Use all physical cores only (with hyper-threading disabled)

## 1. Utilize All Available vCPUs
To run across all 16 vCPUs, update the `Allrun` script by replacing each `runParallel` instance with:

```bash
mpirun -np 16 --use-hwthread-cpus <command> -parallel
```

Don't forget to add the `-parallel` flag after the OpenFOAM command.

This achieves near 100% CPU utilization:

[CPU Usage](../_static/foundation_16_vcpus.png)

## 2. Utilize All Physical Cores
Alternatively, divide your domain into the same number of physical cores and run the simulation.
- The computational resource provides 16 vCPUs, but only 8 physical cores.
- Running with 8 processes (one per core) results in approximately **50%** CPU utilization:

![CPU Usage](../_static/quick-start/system_metrics_50_2tpc.png)

To improve this, you can configure the computational resource so that the number of vCPUs matches the number of physical cores by setting `threads_per_core=1`. Learn more about this setting [here](https://inductiva.ai/guides/how-it-works/machines/hyperthreading). 

The `c2d-highcpu-16` machine has 16 vCPUs by default (`threads_per_core=2`). Setting `threads_per_core=1` reduces the vCPUs to 8 — one per physical core. As a result, running with 8 partitions will fully utilize the CPU, showing **100%** utilization:

![CPU Usage](../_static/quick-start/system_metrics_100.png)

## Results
To evaluate how different configurations impact simulation time and cost, we ran the same case under several machine setups, varying the number of MPI processes and the threading mode (`threads_per_core`). These tests help illustrate how factors like hyper-threading and domain decomposition affect performance.

The table below summarizes the results:

| Machine Type   | Threads per Core | vCPUs Available| MPI Procs |Execution Time | Estimated Cost (USD) |
| -------------- | ---------------- | ---------------|---------- |-------------- | ---------- |
| c2d-highcpu-16 | 2                | 16             |  6        | 2 min, 19s    | 0.0030     |
| c2d-highcpu-16 | 1                | 8              |  8        | 1 min, 58s    | 0.0026     |
| c2d-highcpu-16 | 2                | 16             |  8        | 1 min, 56s    | 0.0025     |
| c2d-highcpu-16 | 2                | 16             |  16       | 1 min, 53s    | 0.0025     |

**Observations**:
- Changing `threads_per_core` from 2 to 1 (i.e., disabling hyper-threading) has minimal impact on performance.Differences are likely within the margin of variation.
- Running with all 16 vCPUs gives the best execution time, though the gains are marginal. Keep in mind that these results reflect a small-scale test case. For larger simulations with more MPI partitions, memory bandwidth and communication overhead can become limiting factors.

> **Note**: Simulation performance is case-dependent. The optimal setup here may not apply to more complex workloads.

Next, we’ll apply the same approach to a larger, more realistic simulation.

## Large-Scale Example: Steady-State CFD Simulation of Wind Flow in the Perdigão Region
To explore how OpenFOAM scales in a more realistic scenario, we ran a steady-state CFD simulation of the [Perdigão region in Portugal](https://journals.ametsoc.org/view/journals/bams/100/5/bams-d-17-0227.1.xml). This area is known for its two parallel ridgelines, which create intricate wind flow patterns and have made it a reference location for atmospheric research.

The simulation was carried out with OpenFOAM’s `simpleFoam` solver on a structured, terrain-following graded mesh containing **14 million cells**. The computational domain covered **30 × 30 × 3 km**, with idealized atmospheric boundary layer (ABL) conditions at the inlet. Turbulence was modeled using the **k–ε closure model**, and we applied a stepped **first- to second-order convection scheme** to improve solution accuracy and convergence.

### Performance Comparison
We ran these tests to understand how hyperthreading and machine size affect OpenFOAM performance and cost. Specifically, we wanted to see whether it was better to:
- (1st row) stick to physical cores on a large machine with hyperthreading enabled,
- (2nd row) try to exploit all vCPUs (logical cores) on the same large machine,
- (3rd row) disable hyperthreading and use only physical cores on the large machine, or
- (4th row) use a smaller machine with fewer physical cores to reduce cost at the expense of runtime.

Below are the results:

| Machine Type   | Threads per Core | vCPUs Available | MPI Procs | Execution Time | Estimated Cost (USD) |
| -------------- | ---------------- | --------------- | --------- | -------------- | ---------- |
| c4d-highcpu-96 | 2                | 96              | 48        | 9h, 20 min   | 15.48      |
| c4d-highcpu-96 | 2                | 96              | 96        | 10h, 58 min  | 18.21      |
| c4d-highcpu-96 | 1                | 48              | 48        | 9h, 23 min   | 15.58      |
| c4d-highcpu-48 | 2                | 48              | 48        | 19h, 8 min   | 15.93      |

### Insights
The first configuration (48 MPI ranks on 96 vCPUs with hyperthreading enabled) follows the standard OpenFOAM-Foundation practice of mapping one rank per physical core. This resulted in the fastest and most cost-effective runtime of 9h 20min and US$15.48.

In the second run, we increased the number of MPI ranks to 96 to fully utilize all available vCPUs. The idea was to test whether fully loading the hyperthreaded machine would improve throughput. Instead, runtime increased to **10h, 58 min**, demonstrating that oversubscribing hyperthreaded cores adds contention and communication overhead, which reduces efficiency.

For the third run, we disabled hyperthreading (`threads_per_core=1`), limiting the machine to 48 vCPUs corresponding exactly to the 48 physical cores. The runtime (**9h, 23 min**) was nearly identical to the first case, confirming that hyperthreading does not provide a performance benefit when already running one MPI rank per physical core.

Finally, we tested a smaller instance (`c4d-highcpu-48`) with 48 vCPUs, which map to only 24 physical cores. Running 48 MPI ranks on this machine oversubscribed the hardware. Although the hourly cost was lower, the runtime increased to 19h 8min, resulting in a slightly higher total cost compared to the larger instance, making this option both slower and more expensive.

**In summary**: the best strategy seems to be running one MPI rank per physical core. Hyperthreading does not provide a measurable benefit for this workload, while trying to use all logical cores or relying on smaller hyperthreaded instances leads to poorer performance and, in some cases, higher costs.