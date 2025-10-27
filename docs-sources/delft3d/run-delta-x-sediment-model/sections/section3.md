# Evaluate the Impact of Hyper-threading
All the simulation runs discussed so far were executed with **hyper-threading enabled**, meaning computational resources were configured with `threads_per_core=2`. This is the **default** setting for virtual machines on Inductiva (learn more [here](https://inductiva.ai/guides/how-it-works/machines/hyperthreading)). 

In traditional HPC environments, however, it's common practice to run one thread per physical core, with **hyper-threading disabled**. This approach helps avoid resource contention and can lead to more predictable and consistent performance.

To disable hyper-threading and ensure only physical cores are used, simply configure your `MachineGroup` with `threads_per_core=1`:

```python
cloud_machine = inductiva.resources.MachineGroup( \
	provider="GCP",
	machine_type="c2d-highcpu-32",
	threads_per_core=1,
	spot=True)
```

The number of available vCPUs will be halved, but the underlying number of physical cores remain the same. 

## Performance with Hyper-threading Disabled
We repeated the **2021 Spring deployment** simulation on the same machine types used before, but this time with hyper-threading disabled. Here are the execution times and associated costs:

| Machine Type      | Hyper-threading  | vCPUs  | Physical Cores   | Execution Time | Estimated Cost (USD) |
|-------------------|------------------|--------|------------------|----------------|----------------------|
| c2d-highcpu-32    | Enabled          | 32     | 16               | 17h, 54 min    | 2.78                 |
| c2d-highcpu-32    | Disabled         | 16     | 16               | 23h, 31 min    | 3.65                 |
| c2d-highcpu-56    | Enabled          | 56     | 28               | 13h, 14 min    | 3.52                 |
| c2d-highcpu-56    | Disabled         | 28     | 28               | 12h, 5 min     | 3.21                 |
| c2d-highcpu-112   | Enabled          | 112    | 56               | 10h, 33 min    | 5.52                 |
| c2d-highcpu-112   | Disabled         | 56     | 56               | 8h, 20 min     | 4.36                 |
| c4d-highcpu-48    | Enabled          | 48     | 24               | 9h, 30 min     | 7.93                 |
| c4d-highcpu-48    | Disabled         | 24     | 24               | 8h, 26 min     | 7.04                 |
| c4d-highcpu-96    | Enabled          | 96     | 48               | 6h, 49 min     | 11.32                |
| c4d-highcpu-96    | Disabled         | 48     | 48               | 5h, 9 min      | 8.55                 |

## Comparing Results: With vs. Without Hyperthreading
Disabling hyper-threading generally improved performance for the Delft3D simulations, particularly on higher-core machines. Across all tested configurations, runs without hyper-threading were consistently faster and slightly more cost-efficient, likely due to reduced contention for shared CPU resources. However, the gains varied depending on the machine type: smaller configurations showed modest improvements, while larger ones benefited more substantially.

Overall, if your workload is heavily compute-bound and scales efficiently across cores, running without hyper-threading can deliver more stable and efficient performance. Still, hyper-threading remains useful for lighter or mixed workloads where I/O or memory operations dominate.

With **Inductiva**, you can easily choose whether to run machines with or without hyper-threading, allowing you to optimize performance for your specific modeling workload.