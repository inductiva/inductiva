# Welcome to the Inductiva API Tutorials

Here you will find in-depth guides explaining the inner workings
of the Inductiva API,as well as how to use it for achieving certain
pratical goals. These tutorials add to the information snippets
available in the [official API documentation](https://docs.inductiva.ai/en/latest/)
by providing more detailed step-by-step explanations and instructions.

## Available Tutorials

* [**Introduction to Inductiva API**](https://docs.inductiva.ai/en/latest/intro_to_api/how_it_works.html).
In this tutorial, we will give you a comprehensive overview of how
the API works, explaining key concepts and how different components
play together. If you never used the API before, we recommend you
read this tutorial.

<div align="center">
   <img width="75%"
    src="./_static/infographic-apifunctionality-fullscreen.svg"
    alt="Inductiva API Usage Flow">
</br>
</div>
<br>

* [**Generating Synthetic Data for training Physics-ML models**](generating-synthetic-data/synthetic-data-generation-1.md).
Synthetic data allows us to train Physics-ML models when you don't
have enough observational data, which often happens in many practical
scenarios. In those cases, we can use physical simulators to mimic
real-world dynamics under various simulated conditions, and produce
data that can help the model to learn the specifics of the scenario
while respecting the underlying physical laws. This tutorial series
will walk you through the process of using the Inductiva API for
generating one of such synthetic training dataset to enable your ML
models to learn complex fluid dynamics.

<div style="display: flex; justify-content:center">
<video width=250 loop muted autoplay preload="auto">
<source src="./_static/generating-synthetic-data/dambreak.mp4" type="video/mp4">
</video>
</div>
<br>

* [**Partial Differential Equations -- Finite Differences and Physics-Informed Neural Networks**](pdes/heat-1-an-introduction.md).
A step by step tutorial on numerical methods to solve Partial
Differential Equations, where we use the 2D Heat Equation as our
working example. We start with the Finite-Differences method,
introduce the main concepts behind Physics-Informed Neural Networks,
and conclude with a Generalized Neuro-Solver that can handle more
complex geometries and varying initial/boundary conditions;

<div style="display: flex; justify-content:center">
<video width=250 loop muted autoplay preload="auto">
<source src="./_static/pdes/cover_slow.mp4" type="video/mp4">
</video>
</div>
<br>

```{toctree}
---
caption: Built-in Simulators
maxdepth: 1
hidden: true
---
simulators/overview
simulators/AmrWind
simulators/CaNS
simulators/DualSPHysics
simulators/FDS
simulators/FVCOM
simulators/GROMACS
simulators/NWChem
simulators/OpenFOAM
simulators/OpenFAST
simulators/QuantumEspresso
simulators/Reef3D
simulators/SCHISM
simulators/SPlisHSPlasH
simulators/SWAN
simulators/SWASH
simulators/XBeach
simulators/CustomImage
```

```{toctree}
---
caption: Quick Recipes
maxdepth: 1
hidden: true
---
how_to/run-parallel_simulations
how_to/manage_computational_resources
how_to/set-up-elastic-machine-group
how_to/set-up-mpi-cluster
how_to/manage-remote-storage
how_to/manage_tasks
how_to/manage_and_retrieve_results
how_to/reuse-files
how_to/run-benchmarks

```

```{toctree}
---
caption: Synthetic Data for Physics-ML
maxdepth: 1
hidden: true
---

generating-synthetic-data/synthetic-data-generation-1.md
generating-synthetic-data/synthetic-data-generation-2.md
generating-synthetic-data/synthetic-data-generation-3.md
generating-synthetic-data/synthetic-data-generation-4.md
generating-synthetic-data/synthetic-data-generation-5.md
generating-synthetic-data/synthetic-data-generation-6.md
```

```{toctree}
---
caption: Partial Differential Equations
maxdepth: 1
hidden: true
---
pdes/heat-1-an-introduction.md
pdes/heat-2-finite-differences.md
pdes/heat-3-PINN.md
pdes/heat-4-neurosolver.md
```
