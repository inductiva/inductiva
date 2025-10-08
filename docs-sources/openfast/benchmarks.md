# Selecting the Right Machine for OpenFAST Simulations

## The Curious Case of OpenFAST
While OpenFAST comprises several modular components, the simulation of a **single wind turbine** operates on a **single thread**. While frameworks like FAST.Farm or MPI allow multiple turbines to run in parallel, the simulation of each turbine remains serial.

As a result, increasing the number of virtual CPUs (vCPUs) assigned to a virtual machine does not improve the runtime of a single-turbine simulation. What matters most is the performance of a single physical core (including clock speed, architecture efficiency, and memory bandwidth). 

Interestingly, consumer desktop CPUs often operate at higher clock frequencies than their cloud-based counterparts. Many desktop processors run at 4–5 GHz, with overclocking pushing speeds beyond 5.5 GHz. In contrast, most cloud machines,including those in Inductiva’s infrastructure, deliver sustained clock speeds around 3 GHz or lower.

Because of this, a single OpenFAST simulation frequently runs faster on a high-end desktop than on a cloud instance. However, if you need to run hundreds or thousands of simulations, deploying large numbers of low-cost cloud machines to run them in parallel becomes far more time-efficient than running them sequentially on a desktop. This is where **Inductiva** offers greater flexibility and the challenge shifts to selecting the most suitable cloud machine type for this highly clock-sensitive, single-threaded workload.

## Benchmark Setup
To evaluate performance across cloud machine types, we used the `5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth` test case. 

> **The original input files can be found [here](https://github.com/OpenFAST/r-test/tree/v4.1.0/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN).**

The software versions used for this test are as follows:

| Component              | Version                               |
|------------------------|---------------------------------------|
| OpenMP                 | 201511                                |
| OpenFAST               | v4.1.0                                |
| gcc, g++, gfortran     | 11.4.0 (Ubuntu 11.4.0-1ubuntu1~22.04) |
| kernel                 | 6.6.12-linuxkit                       |

Bear in mind that these a small use case. The results are even more pronounced on a heavier workload.

## Title?

To demonstrate that OpenFAST’s performance doesn’t scale with the number of cores, we ran the same use case on a local machine and on an c2d-highcpu virtual machine with different vCPUs.

Here are the results:

| Machine Type   | vCPUs | Execution Time | Estimated Cost (USD) |
|----------------|-------|----------------|----------------------|
| c2d-highcpu-2  | 2     | 2 min, 47s     | 0.00062              |
| c2d-highcpu-4  | 4     | 2 min, 36s     | 0.00099              |
| c2d-highcpu-8  | 8     | 2 min, 37s     | 0.0018               |
| c2d-highcpu-16 | 16    | 2 min, 33s     | 0.0033               |
| c2d-highcpu-32 | 32    | 2 min, 33s     | 0.0065               |

As you can see, execution time remained almost identical across all machines, regardless of the number of virtual CPUs.

## Choosing the Right VM: Small and Fast Wins
On Inductiva, the most cost-efficient option for single-turbine simulations is to use lightweight virtual machines with 2 vCPUs backed by a single physical core. 

To determine which VM series offers the best runtime for OpenFAST, we evaluated the following compute-optimized machine families (powered by Google Cloud): c2, c2d, c4, and c4d. These represent a progression of Intel and AMD CPU generations and are known for strong per-core performance.

| Machine Type  | Execution Time | Estimated Cost (USD) |
|---------------|-----------------|---------------------|
| c2d-highcpu-2 | 2 min, 47s      | 0.00062             |
| c2d-standard-2| 2 min, 41s      | 0.00069             |
| c4d-highcpu-2 | 1 min, 52s      | 0.0012              |
| c4-highcpu-2  | 2 min, 46s      | 0.0020              |          

The computational resources are configured with `threads_per_core=2` (hyper-threading enabled), which is the **default** setting for virtual machines on Inductiva (learn more [here](https://inductiva.ai/guides/how-it-works/machines/hyperthreading)).

Conclusion

## Should We Disable Hyperthreading for OpenFAST?

Text about it

| Machine Type  | Execution Time | Estimated Cost (USD) |
|---------------|-----------------|---------------------|
| c2d-highcpu-2 | 2 min, 46s      | 0.00062             |
| c2d-standard-2| 2 min, 42s      | 0.00069             |
| c4d-highcpu-2 | 1 min, 43s      | 0.0011              |
| c4-highcpu-2  | 2 min, 44s      | 0.0020              |  

```{banner_small}
:origin: openfast
```





