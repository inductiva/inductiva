# Evaluating the Impact of Hyper-threading
All the simulation runs discussed so far were executed with **hyper-threading enabled**, meaning computational resources were configured with `threads_per_core=2`. This is the **default** setting for virtual machines on Inductiva (learn more [here](https://inductiva.ai/guides/how-it-works/machines/hyperthreading)). 

In traditional HPC environments, however, it's common practice to run one thread per physical core, with **hyper-threading disabled**. This approach helps avoid resource contention and can lead to more predictable and consistent performance.

To disable hyper-threading and ensure only physical cores are used, simply configure your `MachineGroup` with `threads_per_core=1`:

```python
cloud_machine = inductiva.resources.MachineGroup( \
	provider="GCP",
	machine_type="c4d-highcpu-32",
	threads_per_core=1,
	data_disk_gb=20,
	spot=True)
```

The number of available vCPUs will be halved, but the underlying number of physical cores remain the same. 

## Performance with Hyper-threading Disabled
We repeated the `S2_f5/polynya2D_20211007` scenario tests on the same machine types used before, but this time with hyper-threading disabled. Here are the execution times and associated costs:

| Machine Type      | Hyper-threading | vCPUs | Physical Cores | Execution Time | Estimated Cost (US$) |
|-------------------|------------------|--------|------------------|----------------|----------------------|
| c4d-highcpu-32    | Enabled          | 32     | 16               | 2h, 37 min       | 1.45                 |
| c4d-highcpu-32    | Disabled         | 16     | 16               | 2h, 36 min       | 1.45                 |
| c4d-highcpu-48    | Enabled          | 48     | 24               | 1h, 58 min       | 1.64                 |
| c4d-highcpu-48    | Disabled         | 24     | 24               | 1h, 56 min       | 1.61                 |
| c4d-highcpu-64    | Enabled          | 64     | 32               | 1h, 42 min       | 1.89                 |
| c4d-highcpu-64    | Disabled         | 32     | 32               | 1h, 35 min       | 1.74                 |
| c4d-highcpu-96    | Enabled          | 96     | 48               | 1h, 34 min       | 2.60                 |
| c4d-highcpu-96    | Disabled         | 48     | 48               | 1h, 15 min       | 2.07                 |

## Comparing Results: With vs. Without Hyperthreading
From these results, we observe:
- At lower vCPU counts, the difference in performance between hyper-threaded and non-hyper-threaded runs is minimal.
- At higher vCPU counts, however, disabling hyper-threading leads to better performance and lower cost per hour.

This trend suggests that contention for RAM access becomes a limiting factor when too many threads compete for shared hardware resources. 