# Benchmarks
This benchmark report presents a performance and cost comparison across various CPU and GPU configurations, serving as your trusted guide in selecting the right simulation hardware for your computational CaNS projects.

We benchmark a temporal boundary layer case with stable stratification, the same scenario detailed in our [tutorial](https://inductiva.ai/guides/cans/run-temporal-boundary-layer-case). For benchmarking purposes, we executed 1% of the original CaNS simulation.

## Results 📊
Below is a detailed comparison of execution times and costs across different machine types:

| Machine Type    | vCPUs | GPU                | GPU Count | Duration          | Cost      | Speed-up |
| --------------- | ----- | ------------------ | --------- | ----------------- | --------- | -------- |
| c3d-highcpu-90  | 90    | -                  | 0         | 28 minutes 11 sec | 0.49 US$ | 0.28x    |
| c3d-highcpu-180 | 180   | -                  | 0         | 17 minutes 16 sec | 0.64 US$ | 0.46x    |
| c3d-highcpu-360 | 360   | -                  | 0         | 10 minutes 20 sec | 0.89 US$ | 0.78x    |
| a3-highgpu-1    | 26    | NVIDIA H100 (80GB) | 1         | 8 minutes 1 sec   | 0.63 US$ | 1.00x    |
| a3-highgpu-2    | 52    | NVIDIA H100 (80GB) | 2         | 9 minutes 1 sec   | 1.25 US$ | 0.89x    |
| a3-highgpu-4    | 104   | NVIDIA H100 (80GB) | 4         | 6 minutes 5 sec   | 2.05 US$ | 1.32x    |
| a3-highgpu-8    | 208   | NVIDIA H100 (80GB) | 8         | 6 minutes 8 sec   | 3.99 US$ | 1.30x    |

Speed-ups and cost-efficiency gains are expected to be even more significant when running the full simulation.

The data illustrates that scaling up hardware resources does not always translate into linear performance improvements. In some cases, increased overhead, such as inter-GPU communication or suboptimal resource utilization, can lead to even longer runtimes. Careful benchmarking is therefore essential to find the right balance between speed, cost, and efficiency for your specific computational needs.