# The Inductiva Guide to AMR-Wind

Your resource hub for all things AMR-Wind at Inductiva. Whether you're just starting out or an experienced user, you'll find the resources you need to seamlessly run your AMR-Wind simulations on Cloud machines equipped with hundreds of cores and terabytes of disk space.

Inductiva simplifies research by making high-performance computing more accessible and cost-effective. Use the power of the Cloud to **scale your simulations** and **finish your projects sooner**, while keeping your costs in check! 

## About AMR-Wind
[AMR-Wind](https://github.com/Exawind/amr-wind) is a massively parallel, block-structured adaptive-mesh, incompressible flow solver for wind turbine and wind farm simulations. The primary applications for AMR-Wind are: performing large-eddy simulations (LES) of atmospheric boundary layer (ABL) flows and simulating wind farm turbine-wake interactions using actuator disk or actuator line models for turbines.

## What You'll Find Here
- **Tutorials:** Step-by-step guides to help you learn how to run AMR-Wind through the Inductiva API. From getting started to advanced tutorials, we have you covered.
- **Benchmarks:** A trusted guide to selecting the right simulation hardware for your needs. These benchmarks, conducted using the Inductiva platform, provide insight into how AMR-Wind performs on different hardware configurations.

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
caption: Tutorials
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
caption: Visualization
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
