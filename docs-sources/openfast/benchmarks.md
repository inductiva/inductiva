# Benchmarks

## The curious case of OpenFAST
OpenFAST is a software that operates on a single thread, limiting its ability to utilize the full processing capacity of modern CPUs, which often support multiple threads in parallel. As a result, the performance of OpenFAST simulations is primarily influenced by the CPU clock frequency.

Interestingly, consumer desktop computers typically have higher clock frequencies compared to high-performance machines used in cloud centers. For example, desktop CPUs commonly operate between 4 and 5 GHz (with potential overclocking beyond 5.5 GHz), while cloud machines generally run at clock speeds around 3 GHz or lower.

For a single OpenFAST simulation, a desktop machine is typically faster due to its higher clock speeds, making it the more efficient choice for this specific use case.

## Benchmarking OpenFAST
The following benchmarks highlight how the execution speed of OpenFAST simulations is primarily determined by CPU clock speed, rather than the number of cores available.

For this purpose, we used the `5MW_OC4Semi_WSt_WavesWN` example, an extension of the reference case from the “Definition of a 5-MW Reference Wind Turbine for Offshore System Development",
which can be found on the OpenFAST GitHub repository and is also referenced in the [Run 50 Simulations in Parallel](https://inductiva.ai/guides/openfast/OpenFAST_advanced) tutorial.

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
| Machine Type  | Number of VCPUs | Execution Time | Estimated Cost |
|---------------|-----------------|----------------|-----------|
| n2-highcpu-2  | 2               |34.1 s          |0.00013 US$|
| n2-highcpu-4  | 4               |35.0 s          |0.00031 US$|
| n2-highcpu-8  | 8               |30.4 s          |0.00034 US$|
| n2-highcpu-16 | 16              |30.9 s          |0.00065 US$|
| n2-highcpu-32 | 32              |30.2 s          |0.0012 US$ |

The execution time remained almost identical across all machines, regardless of the number of virtual CPUs.






