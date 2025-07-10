# Benchmarks
This benchmark report presents a performance comparison across various GPU configurations, serving as your trusted 
guide in selecting the right simulation hardware for your computational CaNS projects.

We benchmark a temporal boundary layer with stable stratification case, following the same scenario detailed in our 
[tutorial](https://inductiva.ai/guides/cans/run-temporal-boundary-layer-case). For benchmarking purposes, the case was run 
using a slightly coarser mesh, reducing the number of grid points by 25% along each spatial direction compared to the 
original CaNS simulation.

## Results
The benchmarks cover a range of cloud machines with different GPUs. The reference setup is the most affordable and 
smallest configuration, featuring 4 virtual CPUs (vCPUs) paired with a single **NVIDIA L4 GPU**. Other configurations that were 
tested include more powerful machines with increased CPU counts and higher-performance GPUs, such as the **NVIDIA A100** and **H100**, 
which allow us to evaluate how scaling hardware resources affects simulation speed.

Below is a detailed comparison of execution times and speed-ups across different machine types:

| Machine Type    | vCPUs | GPU            | GPU Count | Execution Time| Estimated Cost (USD) | Speed-up  |
| --------------- | ----- | ---------------| --------- | ------------- | ---------- | --------- |
| g2-standard-4   | 4     | NVIDIA L4      | 1         | 25h, 3 min    | 6.86       | Reference |
| g2-standard-24  | 24    | NVIDIA L4      | 2         | 15h, 55 min   | 10.75      | 1.57x     |
| a2-highgpu-1    | 12    | NVIDIA A100    | 1         | 4h, 44 min    | 7.38       | 5.29x     |
| a2-highgpu-2    | 24    | NVIDIA A100    | 2         | 2h, 47 min    | 8.85       | 9.00x     |
| a3-highgpu-1    | 26    | NVIDIA H100    | 1         | 2h, 26 min    | 6.52       | 10.29x    |
| a3-highgpu-2    | 52    | NVIDIA H100    | 2         | 1h, 36 min    | 8.64       | 15.65x    |

## Summary
The benchmark results clearly demonstrate the substantial performance gains achievable by leveraging more powerful 
GPU configurations for CaNS simulations. The fastest machine tested, equipped with **NVIDIA H100 GPUs**, achieved up to 
a **15.7× speed-up** compared with the smallest GPU tested, a configuration with 4 vCPUs and a single NVIDIA L4 GPU.

With Inductiva, you're able to seamlessly select the hardware that delivers the performance your simulations demand.