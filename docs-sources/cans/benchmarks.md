# GPU Analysis & Results

Benchmarks conducted by **Inductiva** with technical support from **Pedro Costa (TU Delft)**

*Special thanks to **Dr. Baptiste Hardy (TU Delft)** for his support in devising this temporal boundary layer setup*

---

This benchmark report presents a performance comparison across various GPU configurations, serving as your trusted 
guide in selecting the right simulation hardware for your computational CaNS projects.

We benchmark a temporal boundary layer with stable stratification case, following the same scenario detailed in our 
[tutorial](run-temporal-boundary-layer-case).

## Results
The benchmarks cover a range of cloud machines with different GPUs. The reference setup is the most affordable and smallest 
configuration, featuring 4 virtual CPUs (vCPUs) paired with a single **NVIDIA L4 GPU**. Other configurations that were tested 
include more powerful machines with increased CPU counts and higher-performance GPUs, such as the **NVIDIA A100** and **H100**, 
which allow us to evaluate how scaling hardware resources affects simulation speed.

Below is a detailed comparison of execution times and speed-ups across different machine types:

| Machine Type    | vCPUs | GPU Type       | GPU Count | Execution Time | Speed-up  | Estimated Cost (USD)  |
|-----------------|-------|----------------|-----------|----------------|-----------|-----------------------|
| g2-standard-4   | 4     | NVIDIA L4      | 1         | 25h, 3 min     | Reference | 6.86                  |
| g2-standard-24  | 24    | NVIDIA L4      | 2         | 15h, 55 min    | 1.57x     | 10.75                 |
| a2-highgpu-1    | 12    | NVIDIA A100    | 1         | 4h, 44 min     | 5.29x     | 7.38                  |
| a2-highgpu-2    | 24    | NVIDIA A100    | 2         | 2h, 47 min     | 9.00x     | 8.85                  |
| a3-highgpu-1    | 26    | NVIDIA H100    | 1         | 2h, 26 min     | 10.29x    | 6.52                  |
| a3-highgpu-2    | 52    | NVIDIA H100    | 2         | 1h, 36 min     | 15.65x    | 8.64                  |

<p align="center"><em>Table 1: Benchmark results on Inductiva</em></p>

To further validate the performance, we calculated the wall-time per time step, per grid cell, per GPU. Each time step in CaNS 
involves 3 RK3 substeps, with a large Poisson equation being solved at each substep. In this simulation mode of CaNS (with 
`is_impdiff_1d = T` set in the configuration file), the expected performance is around 1 nanosecond per grid cell per GPU. 
These estimates are summarized below:

| Machine Type    | GPU Type     | GPU Count | Execution Time | Time per Cell-Step × GPUs (s) |
|-----------------|--------------|-----------|----------------|-------------------------------|
| g2-standard-4   | NVIDIA L4    | 1         | 25h, 3 min     | 1.769e-08                     |
| g2-standard-24  | NVIDIA L4    | 2         | 15h, 55 min    | 2.247e-08                     |
| a2-highgpu-1    | NVIDIA A100  | 1         | 4h, 44 min     | 3.343e-09                     |
| a2-highgpu-2    | NVIDIA A100  | 2         | 2h, 47 min     | 3.929e-09                     |
| a3-highgpu-1    | NVIDIA H100  | 1         | 2h, 26 min     | 1.719e-09                     |
| a3-highgpu-2    | NVIDIA H100  | 2         | 1h, 36 min     | 2.259e-09                     |

<p align="center"><em>Table 2: Estimates calculated by Pedro Simões Costa</em></p>

## Summary
The benchmark results clearly demonstrate the substantial performance benefits of using higher-end GPUs for CaNS simulations. 
The best-performing setup, equipped with two **NVIDIA H100** GPUs, achieved a **15.7× speed-up** over the baseline machine with 
a single NVIDIA L4 GPU, reducing execution time from 25 hours to just 1 hour and 36 minutes.

These gains, however, are subject to the limits of strong scaling: as GPU count increases while the problem size remains fixed, 
each GPU handles a smaller workload, leading to reduced occupancy and less-than-linear scaling. 

Comparing the measured execution times (Table 1) with theoretical estimates derived from wall-time per grid cell per GPU (Table 2), 
we observe strong agreement across all machines. The anticipated performance of approximately 1 nanosecond per cell per GPU 
was consistently achieved, confirming the robustness and reliability of both the Inductiva platform and the CaNS solver.

With Inductiva, you’re able to seamlessly select the hardware that delivers the performance your simulations demand.

```{banner_small}
:origin: cans
```