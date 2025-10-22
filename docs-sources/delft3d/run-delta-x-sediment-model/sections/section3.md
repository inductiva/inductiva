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

| Machine Type      | Hyper-threading | vCPUs | Physical Cores | Execution Time | Estimated Cost (US$) |
|-------------------|------------------|--------|------------------|----------------|----------------------|


## Comparing Results: With vs. Without Hyperthreading
<Conclusion>