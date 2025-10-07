# Choosing the Right Virtual Machine

## The curious case of OpenFAST
While OpenFAST comprises several modular components, the simulation of a **single wind turbine** operates on a **single thread**. While frameworks like FAST.Farm or MPI allow multiple turbines to run in parallel, the simulation of each turbine remains serial.

As a result, increasing the number of virtual CPUs (vCPUs) assigned to a virtual machine does not improve the runtime of a single-turbine simulation. What matters most is the performance of a single physical core (including clock speed, architecture efficiency, and memory bandwidth). 

Interestingly, consumer desktop CPUs often operate at higher clock frequencies than their cloud-based counterparts. Many desktop processors run at 4–5 GHz, with overclocking pushing speeds beyond 5.5 GHz. In contrast, most cloud machines,including those in Inductiva’s infrastructure, deliver sustained clock speeds around 3 GHz or lower.

Because of this, a single OpenFAST simulation frequently runs faster on a high-end desktop than on a cloud instance. However, if you need to run hundreds or thousands of simulations, deploying large numbers of low-cost cloud machines to run them in parallel becomes far more time-efficient than running them sequentially on a desktop. This is where **Inductiva** offers greater flexibility and the challenge shifts to selecting the most suitable cloud machine type for this highly clock-sensitive, single-threaded workload.













## Benchmarking OpenFAST
The following benchmarks highlight how the execution speed of OpenFAST simulations is primarily determined by CPU clock speed, rather than the number of cores available.

For this purpose, we used the `5MW_OC4Semi_WSt_WavesWN` example, an extension of the reference case from the “Definition of a 5-MW Reference Wind Turbine for Offshore System Development",
which can be found on the OpenFAST GitHub repository and is also referenced in the [Run 50 Simulations in Parallel](run-50-simulations-in-parallel/index) tutorial.

### Software Versions
The software versions used for this benchmark are as follows:

| Component              | Version                               |
|------------------------|---------------------------------------|
| OpenMP                 | 201511                                |
| OpenFAST               | v4.0.2                                |
| gcc, g++, gfortran     | 11.4.0 (Ubuntu 11.4.0-1ubuntu1~22.04) |
| kernel                 | 6.6.12-linuxkit                       |


### Performance and Cost Analysis
To demonstrate that OpenFAST’s performance doesn’t scale with the number of cores, we ran the same use case on a local machine and on an n2-highcpu virtual machine with different vCPUs.

Here are the results:
| Machine Type  | vCPUs | Execution Time | Estimated Cost (USD) |
|---------------|-----------------|----------------|-----------|
| n2-highcpu-2  | 2               |34.1 s          |0.00013|
| n2-highcpu-4  | 4               |35.0 s          |0.00031|
| n2-highcpu-8  | 8               |30.4 s          |0.00034|
| n2-highcpu-16 | 16              |30.9 s          |0.00065|
| n2-highcpu-32 | 32              |30.2 s          |0.0012 |

The execution time remained almost identical across all machines, regardless of the number of virtual CPUs.

```{banner_small}
:origin: openfast
```





