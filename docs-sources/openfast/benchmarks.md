# Selecting the Right Machine for OpenFAST Simulations

## Understanding OpenFAST’s Single-Threaded Nature
While OpenFAST comprises several modular components, the simulation of a **single wind turbine** operates on 
a **single thread**. While frameworks like FAST.Farm or MPI allow multiple turbines to run in parallel, the simulation 
of each turbine remains serial.

As a result, increasing the number of virtual CPUs (vCPUs) assigned to a virtual machine does not improve the 
runtime of a single-turbine simulation. What matters most is the performance of a single physical core (including 
clock speed, architecture efficiency, and memory bandwidth). 

Interestingly, consumer desktop CPUs often operate at higher clock frequencies than their cloud-based counterparts. 
Many desktop processors run at 4-5 GHz, with overclocking pushing speeds beyond 5.5 GHz. In contrast, most cloud 
machines,including those in Inductiva's infrastructure, deliver sustained clock speeds around 3 GHz or lower.

Because of this, a single OpenFAST simulation frequently runs faster on a high-end desktop than on a cloud instance. 
However, if you need to run hundreds or thousands of simulations, deploying large numbers of low-cost cloud machines to run them in parallel becomes far more time-efficient than running them sequentially on a desktop. This is where **Inductiva** offers greater flexibility and the challenge shifts to selecting the most suitable cloud machine type for this highly clock-sensitive, single-threaded workload.

## Simulation Setup
To evaluate performance across different virtual machine types, we used the `5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth` test case. 
This is an offshore, fixed-bottom NREL 5-MW turbine simulation, with the majority of the computational expense occurring in the HydroDyn wave-dynamics calculation, as referenced [here](https://github.com/OpenFAST/r-test/tree/v4.1.0/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN).

The physics modules used in this case are: ElastoDyn, InflowWind, AeroDyn15, ServoDyn, HydroDyn and SubDyn.

The simulation runs for 60 seconds of turbine operation using a time step of 0.01 seconds.

> 🔗 **Input files for this test case are available [here](https://github.com/OpenFAST/r-test/tree/v4.1.0/glue-codes/openfast/5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth).**

The software environment used for this test:

| Component              | Version                               |
|------------------------|---------------------------------------|
| OpenFAST               | v4.1.0                                |
| gcc, gfortran          |(Ubuntu 11.4.0-1ubuntu1~22.04.2) 11.4.0|
| kernel                 |         6.8.0-85-generic              |

Note that this is a relatively lightweight test case. On heavier workloads, performance differences between machine types become even more pronounced.

## OpenFAST Performance vs. vCPU Count
To demonstrate that OpenFAST’s performance is not influenced by the number of allocated vCPUs, we ran the same simulation across multiple c2d-highcpu cloud machines on Inductiva, each with a different number of vCPUs.

| Machine Type   | vCPUs | Execution Time | Estimated Cost (USD) |
|----------------|-------|----------------|----------------------|
| c2d-highcpu-2  | 2     | 2 min, 47s     | 0.00062              |
| c2d-highcpu-4  | 4     | 2 min, 36s     | 0.00099              |
| c2d-highcpu-8  | 8     | 2 min, 37s     | 0.0018               |
| c2d-highcpu-16 | 16    | 2 min, 33s     | 0.0033               |
| c2d-highcpu-32 | 32    | 2 min, 33s     | 0.0065               |

As shown, the **execution time remains effectively constant** regardless of the number of vCPUs. This confirms that OpenFAST only utilizes a single thread, and additional virtual cores simply remain idle. Meanwhile, costs increase linearly with vCPU count, making higher-core machines inefficient for single-turbine runs.

## When Less is More: Selecting the Right Machine
To identify the most effective virtual machines for single-turbine OpenFAST simulations, we evaluated four compute-optimized VM families (powered by Google Cloud): `c2`, `c2d`, `c4`, and `c4d`. All tests used machines with 2 vCPUs, isolating the effect of underlying hardware generation.

| Machine Type    | CPU Type             | Execution Time | Estimated Cost (USD) |
|-----------------|----------------------|----------------|---------------------|
| c2d-highcpu-2   | AMD EPYC             | 2 min, 46 sec  | 0.00062             |
| c2d-standard-2  | AMD EPYC             | 2 min, 42 sec  | 0.00069             |
| c4d-highcpu-2   | AMD EPYC             | 1 min, 43 sec  | 0.0011              |
| c4-highcpu-2    | Intel Xeon Scalable  | 2 min, 44 sec  | 0.0020              |

The computational resources are configured with `threads_per_core=2` (hyper-threading enabled), which is the **default** setting for virtual machines on Inductiva (learn more [here](https://inductiva.ai/guides/how-it-works/machines/hyperthreading)).

The best performer in terms of raw speed is the `c4d-highcpu-2`. This speedup is achieved without increasing the vCPU count, demonstrating that newer-generation CPUs can significantly accelerate single-threaded workloads like OpenFAST, purely due to architectural improvements and higher per-core performance.

That said, the `c4d-highcpu-2` also comes with a higher cost. If your goal is to maximize simulation throughput per dollar, the `c2d-highcpu-2` remains the best value option.

## Disabling Hyper-Threading: Does It Matter?
In many traditional HPC environments, hyper-threading is often disabled to avoid resource contention and to ensure predictable performance for highly parallel, compute-intensive workloads.

Here are the performance results for the same machine types with hyper-threading disabled:

| Machine Type    | CPU Type | Execution Time | Estimated Cost (USD) |
|-----------------|----------|----------------|---------------------|
| c2d-highcpu-2   | AMD EPYC      | 2 min, 46 sec  | 0.00062             |
| c2d-standard-2  | AMD EPYC   | 2 min, 42 sec  | 0.00069             |
| c4d-highcpu-2   | AMD EPYC     | 1 min, 43 sec  | 0.0011              |
| c4-highcpu-2    | Intel Xeon Scalable     | 2 min, 44 sec  | 0.0020              |

Performance tests show **negligible differences** in runtime when hyper-threading is disabled compared to when it is enabled.

For OpenFAST, disabling hyper-threading provides minimal or no benefit. In some cases, enabling hyper-threading can even improve performance by better utilizing available CPU resources.

```{banner_small}
:origin: openfast
```





