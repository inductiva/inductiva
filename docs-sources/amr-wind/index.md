# The Inductiva Guide to AMR-Wind
Your resource hub for all things AMR-Wind at Inductiva. Whether you're just starting out or an experienced user, you'll find the resources you need to seamlessly run your AMR-Wind simulations on Cloud machines equipped with hundreds of cores and terabytes of disk space.

Inductiva simplifies research by making high-performance computing more accessible and cost-effective. Use the power of the Cloud to **scale your simulations** and **finish your projects sooner**, while keeping your costs in check! 

<div style="text-align: center; margin: 20px 0;">
  <div style="font-size: 12px; margin-bottom: 6px;">Try our online Python Editor - run AMR-Wind simulations in your browser</div>
  <a href="https://console.inductiva.ai/editor?simulator_name=amr-wind" 
     style="display: inline-block; width: 55%; padding: 16px 24px; font-size: 14px; font-weight: bold; background-color: var(--playground-button); color: black; text-decoration: none; text-align: center; border-radius: 8px;"
     target="_blank">
    Start Simulating Now
  </a>
</div>

## What You'll Find Here

### Tutorials
Step-by-step guides to help you learn how to run AMR-Wind through the Inductiva API. From getting started to advanced tutorials, we have you covered.

* **Getting Started**
    - [Test Your Inductiva Setup](setup-test)
    - [Run Your First Simulation](quick-start)

* **Advanced Tutorials**
    - [Flow Around a Circular Cylinder](run-flow-cylinder-case)
    - [Run AMR-Wind Simulations Across Multiple Machines Using MPI](mpi-cluster-tutorial)
    - [AMR-Wind Post-Processing](yt-for-post-processing)

### Benchmarks
A trusted guide to selecting the right simulation hardware for your needs. These benchmarks, conducted using the Inductiva platform, provide insight into how AMR-Wind performs on different hardware configurations.

```{banner}
:origin: amr_wind
```

```{toctree}
---
caption: " "
maxdepth: 1
hidden: true
---
versions-and-containers
```

```{toctree}
---
caption: 🛠️ Tutorials
maxdepth: 2
hidden: true
---
setup-test
quick-start
Flow Around a Circular Cylinder <run-flow-cylinder-case>
Run AMR-Wind Simulations Across Multiple Machines <mpi-cluster-tutorial>
```

```{toctree}
---
caption: 📊 Visualization
maxdepth: 2
hidden: true
---
yt for Post-processing <yt-for-post-processing>
```

```{toctree}
---
caption: " "
maxdepth: 1
hidden: true
---
faq
```
