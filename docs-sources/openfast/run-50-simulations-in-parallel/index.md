---
og:image: "https://inductiva.ai/builds/openfast/_static/openfast_animation_30_fps.gif"
---

# Run 50 simulations in parallel
If you only need to run a single OpenFAST simulation, then you should run it on your desktop machine: it will be faster there due to the much higher 
clock speeds. However, if you need to run hundreds or thousands of OpenFAST simulations, you can use the Inductiva API to spin up hundreds of very cheap cloud machines 
to run all these simulations in parallel, instead of running them sequentially on your machine! This is hundreds of times faster. And it is super easy 
(and cost-effective) to do with Inductiva.

In this tutorial, we are going to show you how you can
use the Inductiva API to accelerate your OpenFAST projects, by showing how to
run dozens of simulations in parallel.

<p align="center"><img src="../_static/openfast_animation_30_fps.gif" alt="OpenFAST simulation visualization" width="700"></p>

To demonstrate this, we will use the [`5MW_OC4Semi_WSt_WavesWN`](https://github.com/OpenFAST/r-test/tree/v4.0.2/glue-codes/openfast/5MW_OC4Semi_WSt_WavesWN) example, available on the [OpenFAST GitHub repository](https://github.com/openfast).

This example is an extension of the reference case described in ["Definition of a 5-MW Reference Wind Turbine for Offshore
System Development"](https://www.nrel.gov/docs/fy09osti/38060.pdf).

The tutorial is divided into the following sections:
- [Set up the example files](https://inductiva.ai/guides/openfast/run-50-simulations-in-parallel/sections/section1)
- [Run a single simulation using the Inductiva API](https://inductiva.ai/guides/openfast/run-50-simulations-in-parallel/sections/section2)
- [Generalize the use case](https://inductiva.ai/guides/openfast/run-50-simulations-in-parallel/sections/section3)
- [Launch and manage 50 simulations in parallel on the cloud](https://inductiva.ai/guides/openfast/run-50-simulations-in-parallel/sections/section4)
- [Post-process and analyze simulation results](https://inductiva.ai/guides/openfast/run-50-simulations-in-parallel/sections/section5)

```{banner_small}
:origin: openfast
```

```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
sections/section4.md
sections/section5.md
```
