# Scale Up Your Simulations

Scale Up Your Simulations with Ease

We've prepared a collection of helpful guides to support deeper technical understanding 
for advanced users. This section guides you through benchmarking, monitoring, and optimizing 
your workflows to maximize speed and minimize costs, so you can scale your projects confidently 
and effectively.

## 📘 What’s Inside:
🔢 [Parallel Simulations](parallel-simulations/index) – Learn how to run many simulations at 
once using cloud resources, accelerating large experiments or parameter sweeps with minimal setup.

📁 [Projects](projects/index) – Organize your work into projects to group related simulations, 
streamline task tracking, and manage resource usage efficiently.

⏩ [Recipes](recipes/index) – Ready-to-use code snippets to streamline your workflows.

⏩ [Optimize Workflow](optimize-workflow/index) – Tips and hacks to streamline your workflows.

💰 [Save Costs](save-costs/index) – Tips and strategies to reduce credit usage, from managing 
storage and quotas to selecting efficient machines and cleaning up resources.

🧪 [Generate Dataset](generate-dataset/generate-dataset) – Use simulations to generate custom 
datasets at scale, ideal for machine learning or research requiring labeled physical data.

📊 [Benchmarks](benchmark/index) – Learn to benchmark your simulation on different machines to 
find the most cost-effective and performant configuration.


## 💡 Why It’s Useful:
✓ Reduce costs and simulation time by choosing the most efficient machine type for your task.

✓ Benchmark and optimize your workloads before committing to full-scale runs.

✓ Gain performance insights with system metrics like CPU, memory, and disk usage.

✓ Scale confidently knowing your simulations are running on the right resources. 

✓ Support high-performance use cases, from academic research to commercial engineering, 
with minimal trial and error.


Discover how to efficiently run larger and more complex simulations on Inductiva!   


```{banner}
:origin: scale_up
```


```{toctree}
---
caption: Parallel Simulations
maxdepth: 4
hidden: true
---
Overview <parallel-simulations/index>
Parallel Simulations <parallel-simulations/run-parallel-simulations>
Templating <parallel-simulations/templating>
Parallel Simulations with Templating <parallel-simulations/run-parallel-simulations-with-templating>
Set up an Elastic Machine Group <parallel-simulations/set-up-elastic-machine-group>
Set up an MPI Cluster <parallel-simulations/set-up-mpi-cluster>

```

```{toctree}
---
caption: Projects
maxdepth: 3
hidden: true
---

Overview <projects/index>
projects/projects
projects/manage-projects
projects/visualize-projects

```

```{toctree}
---
caption: Recipes
maxdepth: 3
hidden: true
---
Overview <recipes/index>
♻️ Reuse Files Across Multiple Simulations <recipes/reuse-files>
⬇️ Download specific files from a group of tasks <recipes/download-file-from-project>
🗑️ Clean Up Storage by Condition <recipes/storage-related/index>
👀 Real-Time Monitoring & Conditional Auto Termination <recipes/real-time-simulation-monitoring>
⏰ Setting a Time-to-Live on Your Simulations <recipes/set-task-ttl/set-task-ttl>
```

```{toctree}
---
caption: Optimize Workflow
maxdepth: 2
hidden: true
---
Overview <optimize-workflow/index>
Alerts & Events <optimize-workflow/alerts-events/sections/alerts>

```
 
```{toctree}
---
caption: Save Costs
maxdepth: 2
hidden: true
---
Overview <save-costs/index>
Minimize simulation data <save-costs/save_storage>

```
 
```{toctree}
---
caption: Generate Dataset
maxdepth: 3
hidden: true
---
🧪 Generate a Dataset <generate-dataset/generate-dataset>
```

```{toctree}
---
caption: Benchmarks
maxdepth: 2
hidden: true
---
Overview <benchmark/index>
benchmark/benchmarking
benchmark/why-benchmarks
Run a Benchmark <benchmark/run-benchmarks>
benchmark/monitor-live

```
