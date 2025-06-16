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
    <td>Duration (min:s)</td>
    <td>Speedup</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>1</td>
    <td>112</td>
    <td>144:00</td>
    <td>Baseline</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>2</td>
    <td>224</td>
    <td>90:00</td>
    <td>1.60x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>4</td>
    <td>448</td>
    <td>59:10</td>
    <td>2.43x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>8</td>
    <td>896</td>
    <td>42:12</td>
    <td>3.41x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>12</td>
    <td>1344</td>
    <td>34:15</td>
    <td>4.20x</td>
  </tr>
  <tr>
    <td>c2d-highmem-112</td>
    <td>16</td>
    <td>1792</td>
    <td>30:48</td>
    <td>4.68x</td>
  </tr>
</table>

Runtime decreased as the cluster size increased. With **16 machines** (896 vCPUs), the simulation completed 
in **30 minutes and 48 seconds**, achieving a **speedup of 4.68x** compared to the single-machine runtime of 144 minutes.

Although adding more machines yields noticeable speedups, for a problem of this size, the benefits quickly taper off due to diminishing returns as the cluster size increases.

<br>

> 🛠️ Interested in running AMR-Wind simulations across multiple machines using MPI? Check out this [tutorial](https://inductiva.ai/guides/amr-wind/mpi-cluster-tutorial). 
