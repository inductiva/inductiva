# Welcome to the Inductiva API Tutorials

Here you will find in-depth guides explaining the inner workings
of the Inductiva API,as well as how to use it for achieving certain
pratical goals. These tutorials add to the information snippets
available in the [official API documentation](https://docs.inductiva.ai/en/latest/)
by providing more detailed step-by-step explanations and instructions.

## Available Tutorials

* [**Introduction to Inductiva API**](intro_to_api/how_it_works.md).
In this tutorial, we will give you a comprehensive overview of how
the API works, explaining key concepts and how different components
play together. If you never used the API before, we recommend you
read this tutorial.

<div align="center">
   <img width="75%"
    src="../_static/infographic-apifunctionality-fullscreen.svg"
    alt="Inductiva API Usage Flow">
<\br>
</div>

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
<source src="../_static/generating-synthetic-data/dambreak.mp4" type="video/mp4">
</video>
</div>

```{toctree}
---
caption: Introduction to Inductiva API
maxdepth: 1 
hidden: true
---
intro_to_api/how_it_works
intro_to_api/tasks
intro_to_api/shared_dedicated_resources
intro_to_api/data_flow
intro_to_api/computational-infrastructure
intro_to_api/templating
intro_to_api/configuring-simulators

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
```
