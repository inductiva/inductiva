# CP2K Benchmark: H2O-512 Scaling on the Cloud

<p align="right">
    <small>Last Updated on: 2025-09-09</small>
</p>

## 📌 Overview

This benchmark evaluates the performance of the **CP2K** simulator by running an **ab-initio molecular dynamics (AIMD)** simulation of liquid water.  
The test is based on the standard **H2O-64** benchmark from CP2K but extended to **H2O-512** (512 water molecules, 1536 atoms, 4096 electrons) to better leverage **cloud scalability**.

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

## ⚡ Results

We evaluated performance across multiple configurations, comparing execution times and analyzing the **cost vs runtime** trade-off.

### 1. Execution Times Comparison

This plot shows how the total runtime of the simulation varies across different configurations.

```{raw} html
:file: ./_static/time_graph.html
```

---

### 2. Cost vs Time

This plot compares the **execution cost** against the **simulation runtime**, helping identify the **sweet spot** between **speed** and **cost efficiency**.

```{raw} html
:file: ./_static/time_graph.html
```

---

## 🚀 Conclusions

- **H2O-512** is better suited than **H2O-64** for **cloud-based scaling benchmarks**.
- Execution time decreases **non-linearly** with more compute resources.
- Optimal performance depends on balancing **time-to-results** with **cost**.
- These results provide guidance for selecting the right configuration based on **budget** and **deadline constraints**.

---

```{banner_small}
:origin: cp2k_benchmark-512
```