# NWChem Benchmark: C240 6-31G Scaling on the Cloud

<p align="right">
    <small>Last Updated on: 2025-09-18</small>
</p>

## 📌 Overview

This benchmark evaluates the performance of the **NWChem** simulator by running the **Hybrid density functional calculation on the C240 Buckyball** simulation.  
This NWChem simulation benchmarks the Gaussian basis set DFT module by running a
PBE0 calculation in direct mode on the C240 system using the 6-31G\* basis set
(3600 functions) without symmetry.



> **The original input files can be found [here](https://nwchemgit.github.io/c240_631gs.nw).**


### 🔬 Simulation Details

| **Parameter**            | **Value**                           |
|---------------------------|-------------------------------------|
| Number of atoms           | 240 (all carbon)                    |
| Basis set                 | 6-31G* (cartesian)                  |
| Method                    | DFT (PBE0 functional)               |
| Mode                      | Direct                              |
| Iterations                | 4 (specified)                       |


### Software Versions
The software versions used in this benchmark were the following:

| Component              | Version                               |
|------------------------|---------------------------------------|
| NWChem                  | 7.2.3                                |
| MPI                  | OpenMPI v4.1.2              |
| Operating System       |Ubuntu 22.04.5 LTS|
| kernel                 | 6.8.0-65-generic                     |


## Computational Resources Tested

The benchmark was conducted on the following machine families:

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
- Some configurations were excluded due to **insufficient memory**:  
  - `c3d-highcpu-16`  
  - `c2d-highcpu-16`  
  - local test machines
  - Other setups with low memory
- For the **C2 family**, only the `standard` configuration was tested, as `highcpu` variants are not available.

See the results of the benchmark on the following pages:
- [Execution Times Comparison](exec-time)
- [Cost vs Time](cost-v-time)

```{toctree}
:hidden:
exec-time
cost-v-time
```