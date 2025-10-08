# AMR-Wind Benchmark: Stable Atmospheric Boundary Layer Scaling on the Cloud

<p align="right">
    <small>Last Updated on: 2025-09-18</small>
</p>

## 📌 Overview

This benchmark evaluates the performance of the **AMR-Wind** simulator by running a **stable atmospheric boundary layer** simulation.  
The simulation uses an uniform mesh resolution and directly follows well-established reference studies in the literature.


> **The original input files can be found [here](https://github.com/Exawind/exawind-benchmarks/tree/main/amr-wind/atmospheric_boundary_layer/stable/input_files).**


### 🔬 Simulation Details

| **Category**        | **Parameter**             | **Value**                  |
|----------------------|---------------------------|----------------------------|
| **Domain**           | Size (x,y,z)              | 400 × 400 × 400 m³         |
|                      | Grid resolution           | 512 × 512 × 512            |
| **Time Control**     | Max simulated time        | 36000 s                    |
|                      | Fixed timestep            | 0.0625 s                   |
|                      | CFL factor                | 0.95                       |

### Software Versions
The software versions used in this benchmark were the following:

| Component              | Version                               |
|------------------------|---------------------------------------|
| AMR-Wind                  | 3.7.0                                |
| MPI                  | OpenMPI v4.1.2              |
| Operating System       |Ubuntu 22.04.5 LTS|
| kernel                 | 6.8.0-83-generic                     |


## Computational Resources Tested

The benchmark was conducted on the following machine series:

### Intel-based
- **C2**  
  Intel Xeon Scalable (Cascade Lake)

- **C3**  
  Intel Xeon Scalable (Sapphire Rapids)

- **C4**  
  Intel Xeon Scalable (Granite Rapids) — latest generation

### AMD-based
- **C2D**  
  AMD EPYC (Rome)

- **C3D**  
  AMD EPYC (Genoa)

- **C4D**  
  AMD EPYC (Turin) — latest generation

### Notes on Machine Selection

- This simulation requires approximately **150 GB of memory** (with actual usage varying depending on the number of cores).  
  As a result, some machine configurations were excluded due to **insufficient resources**, including:
  - `c3d-highcpu-60`  
  - `c4d-highcpu-64`  
  - Local test machines  
  - Other similar or lower end setups

- Within the **C2 series**, only the `standard` configuration was evaluated, since `highcpu` variants are not available.


See the results of the benchmark on the following pages:
- [Execution Times Comparison](exec-time)
- [Cost vs Time](cost-v-time)

```{toctree}
:hidden:
exec-time
cost-v-time
```