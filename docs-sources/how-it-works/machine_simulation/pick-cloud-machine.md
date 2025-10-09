# Pick the Right Cloud Machine

Choosing the right machine to run your simulations can make the difference between a simulation that runs in hours versus days, and between spending $10 versus $100 on compute costs. With hundreds of machine types available through Inductiva, this guide will help you navigate the options and make informed decisions that optimize both performance and cost.

## Quick Start

- **Browse all options:** View our [complete machine catalog](https://console.inductiva.ai/machine-groups/instance-types) to browse all available machine type options.

- **Make data-driven decisions:** Run [benchmarks](https://inductiva.ai/guides/scale-up/benchmark/index) on candidate machines with representative workloads.

- **Optimize costs:** Use [spot machines](../machines/spot-machines.md) for up to 60% savings.

## Machine series

Inductiva provides access to Google Cloud Platform machine series, each optimized for specific computational patterns:

### Compute-Optimized Machines
**Best for:** CPU-intensive simulations, mathematical modeling, fluid dynamics
- **series:** C2, C2D, H3
- **Key strength:** Highest performance per core with premium processors
- **Memory ratio:** Optimized compute-to-memory ratios (typically 2-4GB RAM per vCPU)

[Learn more about compute-optimized machines →](https://cloud.google.com/compute/docs/compute-optimized-machines)

### Memory-Optimized Machines
**Best for:** Large-scale simulations requiring extensive data in memory
- **series:** M3
- **Key strength:** High memory-to-vCPU ratios (8GB+ per vCPU)
- **When to choose:** Your simulation loads large datasets or requires extensive intermediate storage

[Learn more about memory-optimized machines →](https://cloud.google.com/compute/docs/memory-optimized-machines)

### General-Purpose Machines
**Best for:** Balanced workloads, development, small to medium simulations
- **series:** C3, C3D, C4,
- **Key strength:** Versatile balance of compute, memory, and networking
- **Typical use cases:** Prototyping, mixed workloads, cost-sensitive applications
- **Cost advantage:** Best price-performance ratio for varied workloads

[Learn more about general-purpose machines →](https://cloud.google.com/compute/docs/general-purpose-machines)

### Accelerator-Optimized Machines
**Best for:** GPU-accelerated simulations, machine learning, parallel computing
- **series:** A2, A3, G2
- **Key strength:** High-performance GPUs (NVIDIA Tesla, A100, H100)
- **Performance boost:** 10-100x speedup for compatible workloads

[Learn more about accelerator-optimized machines →](https://cloud.google.com/compute/docs/accelerator-optimized-machines)

## Machine Naming Convention

Machine names follow a consistent pattern: `series-profile-vcpus`

### Series Identifiers
- **Generation:** C**4** (2024), C**3** (2023) A**2**(2020)
- **Architecture:**
  - **D suffix** → AMD processors (C2D, N2D)
  - **A suffix** → ARM processors (C4A, T2A)
  - **G/A prefix** GPU-equipped (G2, A3)
  - **Default** Intel processors (C2, N2, C3)

### Memory Profiles
- **highcpu:** 1-3GB RAM per vCPU
- **standard:** 3-7GB RAM per vCPU
- **highmem:** 7-14GB RAM per vCPU
- **megamem:** 14-19GB RAM per vCP
- **ultramem:** 24-31GB RAM per vCPU

**Example** `c3d-standard-96` = C3 series, AMD processors, standard memory, 96 vCPUs

## Workload Matching

### CPU-Intensive Workloads
Mathematical computations, iterative algorithms, single-threaded bottlenecks.

**Examples:** Finite element analysis, ray tracing, numerical modeling  
**Recommended:** Compute-optimized series (C2, H3)  
**Memory profile:** `highcpu` or `standard`

### Memory-Intensive Workloads
Large datasets, extensive meshes, in-memory processing requirements.

**Examples:** Large CFD simulations, molecular dynamics, structural analysis  
**Indicator:** >8GB RAM per CPU core requirement  
**Recommended:** Memory-optimized series or `highmem` variants

### Parallel Workloads
MPI-based applications, embarrassingly parallel problems, multi-threaded code.

**Examples:** Monte Carlo simulations, parameter sweeps, parallel finite element  
**Recommended:** High vCPU count machines (96+ cores)  
**Consider:** Multiple smaller machines vs. single large machine

### GPU-Accelerated Workloads
CUDA, OpenCL, or specialized parallel processing frameworks.

**Examples:** Machine learning simulations, GROMACS, custom CUDA applications  
**Recommended:** A2, A3, or G2 series  
**Key factor:** GPU memory capacity matching problem size

## Performance Comparison

Let's compare the specs of the machines we make available with a familiar case, such as a laptop with 32GB of RAM and a 16 core CPU, similar to one that you may own.

At a first glance, such a laptop seems like a reasonably powerful machine and you might be using a similar machine to run relatively large simulations, although sometimes you may have to leave it crunching number for a couple of days.

### From your laptop to the cloud

**c3d-standard-30** provides 30 vCPUs and 120GB RAM - nearly doubling your processing power and quadrupling your memory. This represents a conservative first step into cloud computing while delivering meaningful performance improvements.

**c3d-highmem-96** scales up to 96 vCPUs with 768GB RAM - six times more cores and 24 times more memory than your laptop. This configuration can handle significantly larger problems and should deliver substantial performance improvements for most workloads.

**c4-standard-192** represents cutting-edge performance with 192 vCPUs and 768GB RAM using the latest processor generation. With 12 times more cores than your laptop, this machine can tackle computationally intensive problems that would be impractical on desktop hardware.

**c3d-highmem-360** provides massive computational resources with 360 vCPUs and 2880GB RAM - that's 90 times more memory than your laptop. Very likely, this machine will allow you to run your simulations **5 to 10x faster** than in your laptop, and it will, for sure, let you tackle much larger simulation.

> The actual performance improvement you'll see depends heavily on your specific simulation code, how well it parallelizes, and whether your workload is CPU-bound, memory-bound, or I/O-bound. We recommend running [benchmarks](https://inductiva.ai/guides/scale-up/benchmark/index) to measure real performance with your specific workloads.

> When in doubt, use the c2d-highcpu-112 machine: it has 112 vCPUs and 224GB of RAM, providing an impressive cost/performance ratio (price around 4.62 US\$/hour and spot price of just ~0.68 US\$/hour).

```{banner_small}
:origin: pick-cloud-machine
```