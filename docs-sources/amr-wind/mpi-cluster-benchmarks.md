
## Results

Here are the results of running this simulation on multiple MPI cluster compared
with the baseline of running in a single node configuration.

| Machine Type    | Num Machines | vCPUs | Duration (s) | Speedup |
| --------------- | ------------ | ----- | ------------ | ------- |
| c2d-highcpu-112 | 1            | 112   | 1472         | 1.00x   |
| c2d-highcpu-112 | 2            | 224   | 1404         | 1.05x   |
| c2d-highcpu-112 | 4            | 448   | 1461         | 1.01x   |
| c2d-highcpu-112 | 8            | 896   | 977          | 1.51x   |
| c2d-highcpu-112 | 12           | 1344  | 1408         | 1.05x   |
| c2d-highcpu-112 | 16           | 1792  | 1001         | 1.47x   |


For this simulation setup, the **fastest configuration** was achieved using a
cluster with **8 machines (896 vCPUs)**, completing in **977 seconds**.
Interestingly, increasing the number of machines beyond this point did
**not lead to better performance**. This can be attributed to the diminishing
returns of parallelization for smaller or moderately sized problems, where
communication overhead between nodes begins to offset the gains from additional
compute resources.

## What Happens When We Increase the Simulation Complexity?

To investigate how performance scales with problem size, we increased the
simulation complexity by refining the grid resolution. Specifically, we changed:

```
amr.n_cell              = 512 512 184 # Grid cells at coarsest AMRlevel
```

to

```
amr.n_cell              = 1024 1024 368 # Grid cells at coarsest AMRlevel
```

This effectively **increased the total number of grid cells by a factor of 8**,
leading to a significantly more computationally demanding simulation. In these
cases, **larger clusters with more vCPUs are expected to perform better**, as
the computational workload can be more effectively distributed across machines.
The communication overhead becomes relatively smaller compared to the total
compute time, which improves scalability.

We'll now compare the performance across the same cluster configurations under
this higher-resolution setup.

| Machine Type    | Num Machines | Cores | Duration (s) | Speedup |
| --------------- | ------------ | ----- | ------------ | ------- |
| c2d-highmem-112 | 1            | 112   | 8640         | 1.00x   |
| c2d-highmem-112 | 2            | 224   | 5400         | 1.60x   |
| c2d-highmem-112 | 4            | 448   | 3550         | 2.43x   |
| c2d-highmem-112 | 8            | 896   | 2532         | 3.41x   |
| c2d-highmem-112 | 12           | 1344  | 2055         | 4.20x   |
| c2d-highmem-112 | 16           | 1792  | 1848         | 4.68x   |


The results show that, with the higher-resolution grid, simulation runtimes
decreased consistently as more machines were added, unlike in the lower-resolution
case, where gains plateaued beyond 8 machines. The best performance was observed
using 16 machines (1792 vCPUs), bringing the runtime down to just 1848 seconds,
a significant improvement compared to the single-machine setup at 8640 seconds.

However, the speedup does not scale linearly with the number of vCPUs. While the
4 and 8 machine setups provided solid performance gains, adding more machines
beyond that point yielded diminishing returns.
