# Performance Comparison: CPU vs. GPU Configurations
This benchmark report presents a performance and cost comparison across various CPU and GPU configurations. We identify the most effective hardware setup for running a computational fluid dynamics (CFD) workload, balancing speed, cost-efficiency, and resource utilization.

We benchmark the *flow around a circular cylinder case* from the [AMR-Wind GitHub repository](https://github.com/Exawind/amr-wind/tree/main/test/test_files/ib_cylinder_Re_300), simulating flow over a cylinder at a Reynolds number of 10,000 using a 20-million-cell mesh.

> 🛠️ Learn how to run AMR-Wind simulations on Inductiva in this [tutorial](https://inductiva.ai/guides/amr-wind/quick-start). 

# Benchmark Results 📊
Below is a comparison of execution time and cost across different machine types:

<table>
    <tr>
        <td>Platform</td>
        <td>Machine Type</td>
        <td>Machines</td>
        <td>Cores</td>
        <td>GPUS</td>
        <td>Duration</td>
        <td>Cost (USD)</td>
    </tr>
    <tr>
        <td rowspan="4">AMD EPYC Milan</td>
        <td>c2d-highcpu-16</td>
        <td>1</td>
        <td>16</td>
        <td>0</td>
        <td>2 hours, 8 minutes</td>
        <td>0.23 US$</td>
    </tr>
    <tr>
        <td>c2d-highcpu-32</td>
        <td>1</td>
        <td>32</td>
        <td>0</td>
        <td>1 hour, 14 minutes</td>
        <td>0.26 US$</td>
    </tr>
    <tr>
        <td>c2d-highcpu-56</td>
        <td>1</td>
        <td>56</td>
        <td>0</td>
        <td>53 minutes, 29 seconds</td>
        <td>0.32 US$</td>
    </tr>
    <tr>
        <td>c2d-highcpu-112</td>
        <td>1</td>
        <td>112</td>
        <td>0</td>
        <td>27 minutes, 48 seconds</td>
        <td>0.33 US$</td>
    </tr>
    <tr>
        <td rowspan="4">NVIDIA L4 (7680 CUDA, 240 Tensor Cores)</td>
        <td>g2-standard-4</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <td>55 minutes, 26 seconds</td>
        <td>0.26 US$</td>
    </tr>
    <tr>
        <td>g2-standard-24</td>
        <td>1</td>
        <td>2</td>
        <td>2</td>
        <td>40 minutes, 43 seconds</td>
        <td>0.54 US$</td>
    </tr>
    <tr>
        <td>g2-standard-48</td>
        <td>1</td>
        <td>4</td>
        <td>4</td>
        <td>25 minutes, 7 seconds</td>
        <td>0.66 US$</td>
    </tr>
    <tr>
        <td>g2-standard-96</td>
        <td>1</td>
        <td>8</td>
        <td>8</td>
        <td>18 minutes, 14 seconds</td>
        <td>0.96 US$</td>
    </tr>
<tr>
        <td rowspan="4">NVIDIA A100 (6912 CUDA, 432 Tensor Cores)</td>
        <td>a2-highgpu-1g</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <td>20 minutes, 42 seconds</td>
        <td>0.51 US$</td>
    </tr>
    <tr>
        <td>a2-highgpu-2g</td>
        <td>1</td>
        <td>2</td>
        <td>2</td>
        <td>24 minutes, 50 seconds</td>
        <td>1.23 US$</td>
    </tr>
    <tr>
        <td>a2-highgpu-4g</td>
        <td>1</td>
        <td>4</td>
        <td>4</td>
        <td>17 minutes, 56 seconds</td>
        <td>1.78 US$</td>
    </tr>
    <tr>
        <td>a2-highgpu-8g</td>
        <td>1</td>
        <td>8</td>
        <td>8</td>
        <td>16 minutes, 45 seconds</td>
        <td>3.33 US$</td>
    </tr>
    <tr>
        <td rowspan="2">NVIDIA H100 (14,592 CUDA Cores)</td>
        <td>a3-highgpu-1g</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <td>13 minutes, 59 seconds</td>
        <td>0.58 US$</td>
    </tr>
    <tr>
        <td>a3-highgpu-2g</td>
        <td>1</td>
        <td>2</td>
        <td>2</td>
        <td>16 minutes, 15 seconds</td>
        <td>1.36 US$</td>
    </tr>
</table>

## How to Choose the Best Machine for Your CFD Workload
The ideal compute instance for you depends on your priorities. 

**If Speed Is Your Top Priority**: Opt for the fastest GPUs like the NVIDIA H100 or A100 instances to minimize simulation runtime. While these come at a higher hourly cost, the overall faster completion can reduce total expense.

**If Cost Efficiency Is Most Important**: High-core-count CPUs such as the **AMD c2d-highcpu-112** provide solid performance for the lowest cost. Consider entry-level GPUs like the **NVIDIA L4** for some acceleration without significant cost increases.

**If You Are Exploring or Just Getting Started**: Begin with modest configurations to profile your workload and use our benchmarks as a reference to make informed scaling decisions or decide whether GPU acceleration is worthwhile.

Whatever your priorities, Inductiva takes care of the infrastructure, allowing you to focus on what matters most: your CFD projects.