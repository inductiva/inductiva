# Evaluating AMR-Wind Performance on Multi-Machine Clusters
This page presents detailed performance benchmarks of AMR-Wind simulations run across various multi-node MPI cluster configurations. We compare runtimes and speedups against a baseline single-node setup, and explore how increasing simulation complexity affects scalability.

We benchmark the *neutral Atmospheric Boundary Layer case* from the [AMR-Wind GitHub repository](https://github.com/Exawind/exawind-benchmarks/tree/main/amr-wind/atmospheric_boundary_layer/neutral/input_files).

> 🛠️ Interested in running AMR-Wind simulations across multiple machines using MPI? Check out this [tutorial](https://inductiva.ai/guides/amr-wind/flow-cylinder). 

## Benchmark Results: Baseline Simulation
The following table summarizes simulation runtimes and speedups for the baseline configuration, run on clusters of varying sizes:

<table>
  <tr>
    <td>Machine Type</td>
    <td>Nº of Machines</td>
    <td>vCPUs</td>
    <td>Duration (s)</td>
    <td>Speedup</td>
  </tr>
  <tr>
    <td>c2d-highcpu-112</td>
    <td>1</td>
    <td>112</td>
    <td>1472</td>
    <td>1.00x</td>
  </tr>
  <tr>
    <td>c2d-highcpu-112</td>
    <td>2</td>
    <td>224</td>
    <td>1404</td>
    <td>1.05x</td>
  </tr>
  <tr>
    <td>c2d-highcpu-112</td>
    <td>4</td>
    <td>448</td>
    <td>1461</td>
    <td>1.01x</td>
  </tr>
  <tr>
    <td>c2d-highcpu-112</td>
    <td>8</td>
    <td>896</td>
    <td>977</td>
    <td>1.51x</td>
  </tr>
  <tr>
    <td>c2d-highcpu-112</td>
    <td>12</td>
    <td>1344</td>
    <td>1408</td>
    <td>1.05x</td>
  </tr>
  <tr>
    <td>c2d-highcpu-112</td>
    <td>16</td>
    <td>1792</td>
    <td>1001</td>
    <td>1.47x</td>
  </tr>
</table>

The fastest runtime (977 seconds) was achieved with an **8-machine cluster (896 vCPUs)**. Increasing beyond 8 machines did not further improve performance, highlighting the overhead costs of communication between nodes in moderate-sized simulations.

## Scaling with Increased Simulation Complexity
To evaluate scalability under heavier workloads, we increased the grid resolution from:

```
amr.n_cell              = 512 512 184 # Grid cells at coarsest AMRlevel
```

to:

```
amr.n_cell              = 1024 1024 368 # Grid cells at coarsest AMRlevel
```

This change increases the **total number of grid cells** by a **factor of 8**, significantly boosting computational demands. Under this higher-resolution setup, we expect larger clusters with more vCPUs to better distribute the workload, improving scalability as communication overhead becomes comparatively less impactful.

The performance results are summarized below:

<table>
  <tr>
    <td>Machine Type</td>
    <td>Nº of Machines</td>
    <td>Cores</td>
    <td>Duration (s)</td>
    <td>Speedup</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>1</td>
    <td>112</td>
    <td>8640</td>
    <td>1.00x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>2</td>
    <td>224</td>
    <td>5400</td>
    <td>1.60x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>4</td>
    <td>448</td>
    <td>3550</td>
    <td>2.43x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>8</td>
    <td>896</td>
    <td>2532</td>
    <td>3.41x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>12</td>
    <td>1344</td>
    <td>2055</td>
    <td>4.20x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>16</td>
    <td>1792</td>
    <td>1848</td>
    <td>4.68x</td>
  </tr>
</table>

Unlike the baseline case, runtime consistently decreased as the cluster size grew. The best performance was recorded with 16 machines (1792 vCPUs), completing the simulation in 1848 seconds — a major improvement over the single-machine runtime of 8640 seconds.

However, the speedup does not scale linearly with the number of vCPUs. While the 4 and 8 machine setups showed solid performance gains, adding more machines beyond that point yielded diminishing returns.
