# Synthetic Data for Physics-informed Machine Learning
In this tutorial we’ll show you how to use the Inductiva API to generate synthetic data, at scale and in an economic way, for training Physics-ML models.

We’re going to let you in on our very own recipe for generating synthetic data at scale, giving you an overview of the whole process building on an example from a published study. 

This guide is a resource for both machine learning engineers and enthusiasts alike with a thorough understanding of Physics-ML and simulation software.

Generating synthetic data with the Inductiva API is pretty straightforward and involves just a few steps. These include developing an initial simulation model for the “base case” you wish to study, to running the corresponding simulator across thousands of variations so we can create generate a diverse and large-enough dataset.

Our recipe goes as follows:

Set Up the Base Case: We begin with preparing the configuration files for a “base case” simulation model of the system under study. This step demands a thorough understanding of the system and simulation software, often requiring collaboration with specialists in the field. Machine learning engineers will likely need to partner with domain expertise and experience with specific simulation software. Our API currently supports a growing range of simulators covering diverse simulation needs.
Generalize the Base Case: Then, we need to “generalize” the configuration files of the initial base case in such a way that they can be used to describe a large number of variations of such a base case. Our templating mechanism elegantly replaces templated variables in the configuration files with values at runtime. This means that we can easily tweak the simulation (hyper)parameters through Python scripting, and efficiently cover the configuration space of our simulation, using exhaustive or randomized strategies.
Benchmark Computational Resources: Because running a large number of high-fidelity simulations can be a long and costly process, we need to find a good balance between data resolution, time and costs before generating data at scale. Our API provides all the tools needed for this benchmarking process , allowing us to run variations of our “base case” across different hardware and at different levels of fidelity with only a few lines of code.
Generate Synthetic Data in Bulk: Finally, we use our API to fire up hundreds of cloud machines, running thousands of variations of our “base case”. Once the simulations complete, we collect the output data for post-processing. And just like that, we’re all set!

 let’s take a look at a study done by a group of researchers who used data generated from Smoothed Particle Hydrodynamics (SPH) solvers to train a Graph Neural Network (GNN) model.



```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
sections/section4.md
sections/section5.md
```
