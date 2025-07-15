## Why Benchmark

### Understanding Benchmarks
Benchmarking is essential for making informed decisions about computational resources. When you launch a simulation, its performance and cost depend heavily on how well the available hardware — CPU cores, memory, and GPUs are used. Without systematic testing, you risk paying for resources you don’t use efficiently or waiting longer than necessary for results.

A benchmark lets you test different machine types and configurations to find the sweet spot: the setup that gives you the best performance for the lowest cost. Simply adding more cores doesn't always speed things up linearly. By running a benchmark, you can see exactly how your simulation scales and avoid wasting money on hardware that doesn't provide a real benefit.

### When to Benchmark
You should make benchmarking a regular part of your workflow. Consider it whenever you are:
- Choosing between different machine types
- Launching a new type of simulation.
- Significantly changing your model's size or complexity
- Evaluating the impact of simulation parameters on performance
- Validating that performance improvements justify increased costs

### A Practical Approach
The easiest way to start is by running a shorter, representative version of your full simulation across a few different machine configurations. After this short run, check the performance statistics in the [Web Console](https://console.inductiva.ai/benchmarks) to see how efficiently your simulation used each machine or configuration. Based on these insights, you can confidently choose the optimal configuration for your full-scale workload, knowing you're making the most of your time and budget.

In short, a quick benchmark today saves you from unnecessary spending and delays tomorrow.
