# Evaluating AMR-Wind Performance on Multi-Machine Clusters
This page presents detailed performance benchmark of AMR-Wind simulations conducted on various multi-node MPI cluster configurations. Runtimes and speedups are compared against a baseline single-node setup.

For this purpose, we used the *neutral Atmospheric Boundary Layer case* from the [AMR-Wind GitHub repository](https://github.com/Exawind/exawind-benchmarks/tree/main/amr-wind/atmospheric_boundary_layer/neutral/input_files).

## Results
Below are the results of running this simulation on a multi-node MPI cluster, compared to the single-node configuration (shown in the first row):

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
    <td>Baseline</td>
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
</table>

Runtime consistently decreased as the cluster size grew. The best performance was recorded with **8 machines** (896 vCPUs), completing the simulation in **42 minutes and 48 seconds** — a major improvement over the single-machine runtime of **144 minutes**. 

Runtime consistently decreased as the cluster size increased. With **8 machines** (896 vCPUs), the simulation completed in 30 minutes and 48 seconds, achieving a speedup of 4.68 times compared to the single-machine runtime of 144 minutes.

<br>

> 🛠️ Interested in running AMR-Wind simulations across multiple machines using MPI? Check out this [tutorial](https://inductiva.ai/guides/amr-wind/mpi-cluster-tutorial). 
