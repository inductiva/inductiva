# CP2K Benchmark: H2O-512 Scaling on the Cloud

<p align="right">
    <small>Last Updated on: 2025-09-17</small>
</p>

## 📌 Overview

This benchmark evaluates the performance of the **CP2K** simulator by running an **ab-initio molecular dynamics (AIMD)** simulation of liquid water.  
The test is based on the standard **H2O-64** benchmark from CP2K but extended to **H2O-512** (512 water molecules, 1536 atoms, 4096 electrons) to better leverage **cloud scalability**.


> **The original input files can be found [here](https://github.com/cp2k/cp2k/blob/master/benchmarks/QS/H2O-512.inp).**


### 🔬 Simulation Details

| **Property**                | **H2O-512** Configuration                |
|----------------------------|-----------------------------------------|
| **Number of molecules**    | 512                                    |
| **Total atoms**            | 1536                                   |
| **Total electrons**        | 4096                                   |
| **Basis set**             | TZV2P                                  |
| **Plane-wave cutoff**      | 280 Ry                                |
| **Exchange-correlation**   | LDA                                   |
| **MD steps**              | 10                                     |
| **Initial guess**         | Atomic orbitals                        |

### Why H2O-512?

The **H2O-64** benchmark is designed for local runs or small clusters.  
However, when scaling to **high-performance cloud environments**, its small size underutilizes resources and fails to stress-test parallel performance.  
Thus, we use **H2O-512** to:

- Test **scaling efficiency** on multiple cores
- Measure **cost-performance trade-offs** on the cloud  
- Better simulate real-world production workloads

---

### Software Versions
The software versions used in this benchmark were the following:

| Component              | Version                               |
|------------------------|---------------------------------------|
| cp2k                  | 2025.1                                |
| MPI                  | OpenMPI v4.1.6              |
| Operating System       |Ubuntu 22.04.5 LTS|
| kernel                 | 6.8.0-65-generic                     |


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
- Some configurations were excluded due to **insufficient memory**:  
  - `c3d-highcpu-16`  
  - `c2d-highcpu-16`  
  - local test machines  
- For the **C2 series**, only the `standard` configuration was tested, as `highcpu` variants are not available.

See the results of the benchmark on the following pages:
- [Execution Times Comparison](exec-time)
- [Cost vs Time](cost-v-time)

```{toctree}
:hidden:
exec-time
cost-v-time
```