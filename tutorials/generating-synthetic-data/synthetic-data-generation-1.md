---
myst:
  html_meta:
    description: "Learn how to generate large-scale synthetic datasets for Physics-ML models using the Inductiva API, starting with our unique data generation recipe."
    keywords: "Inductiva API, Programming, HPC, Simulation, Tutorial, Synthetic Data Generation, Physics-ML, SPH" 
---

# Introduction

In this tutorial series we'll show you how to use the Inductiva API to generate
synthetic data, at scale and in an economic way, for training Physics-ML models.

We're going to let you in on our very own recipe for generating synthetic data at scale, giving you an overview of the whole process building on an example from a published study. Over the next steps of this tutorial series, we'll break it all down to show you how it's done. This guide is a resource for both machine learning engineers and enthusiasts alike with a thorough understanding of Physics-ML and simulation software.

## Inductiva's Recipe

Generating synthetic data with the Inductiva API is pretty straightforward and
involves just a few steps. These include developing an initial simulation model for the "base case" you wish to study, to running the corresponding simulator across thousands of variations so we can create generate a diverse and large-enough dataset.

Our recipe goes as follows:

1. **Set Up the Base Case:** We begin with preparing the configuration files for a "base case" simulation model of the system under study. This step demands a thorough understanding of the system and simulation software, often requiring collaboration with specialists in the field. Machine learning engineers will likely need to partner with domain expertise and experience with specific simulation software. Our API currently supports a growing range of [simulators](https://tutorials.inductiva.ai/simulators/overview.html)
covering diverse simulation needs.

2. **Generalize the Base Case:** Then, we need to "generalize" the configuration files of the initial base case in such a way that they can be used to describe a large number of variations of such a base case. Our [templating mechanism](https://tutorials.inductiva.ai/intro_to_api/templating.html) elegantly replaces templated variables in the configuration files with values at runtime. This means that we can easily tweak the simulation (hyper)parameters through Python scripting, and efficiently cover the configuration space of our simulation, using exhaustive or randomized strategies.

3. **Benchmark Computational Resources:** Because running a large number of high-fidelity simulations can be a long and costly process, we need to find a good balance between data resolution, time and costs before generating data at scale. Our API provides all the tools needed for this benchmarking process , allowing us to run variations of our "base case" across different hardware and at different levels of fidelity with only a few lines of code.

4. **Generate Synthetic Data in Bulk:** Finally, we use our API to fire up hundreds of cloud machines, running thousands of variations of our "base
case". Once the simulations complete, we collect the output data for
post-processing. And just like that, we're all set!

We'll get into the details of each of these steps in the coming posts of this
tutorial series. But first, let's take a look at a study done by a group of
researchers who used data generated from Smoothed Particle Hydrodynamics (SPH)
solvers to train a Graph Neural Network (GNN) model.

## Learning Complex Physics: A Practical Study

Let's take a closer look at a practical application of synthetic datasets
through this [study by Sanchez-Gonzalez et
al.](https://arxiv.org/abs/2002.09405). The researchers explore how a Graph
Neural Networks (GNNs) model, trained on synthetic data generated using Smoothed Particle Hydrodynamics (SPH) solvers, can be applied to simulate fluid dynamics. [Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) (SPH), stands out for its mesh-free approach that models fluid flows by allowing coordinates to move along with the fluid. It is widely celebrated for its ability to simulate exceptionally realistic fluid dynamics, such as those we see in cinematic visual effects.

The research abstract points out:

> This work uses a state-of-the-art deep learning model to learn the dynamics of fluids and to simulate how they flow without the need for solving the physics equations that govern the phenomenon. The proposed model is based on a Graph Network architecture, which takes as input a graph whose nodes represent fluid particles and edges physical interactions among them. The model performs sequential steps of computations involving the nodes that are connected by edges to predict how the fluid flows. (...) The model is trained on data generated with Smoothed Particle Hydrodynamics, a Lagrangian method in Computational Fluid Dynamics that simulates fluids as a set of particles.

The researchers behind this study consider a simple simulation case: modeling the behavior of a 3D fluid block, made up of many tiny uniformly dense particles all packed together in a cube shape, suspended in mid-air and
left to fall under the effect of gravity inside a box. As the simulation
progresses, the fluid block collides with the walls of the box, resulting in
complex fluid behaviors and surface patterns.

<div style="display: flex; justify-content:center">
<video width=500 loop muted autoplay preload="auto">
<source src="../_static/generating-synthetic-data/dambreak.mp4" type="video/mp4">
</video>
</div>

Video 1: Dynamics of a water block simulated with 32984 SPH
particles. Here, the block is left under the effect of gravity - *splash!*
Simulation performed via Inductiva API.

The core challenge here lies in training a Physics-ML model to predict the
velocity of each fluid particle over time while adhering to the laws of fluid
dynamics that are not explicitly encoded within the ML model. This
involves generalizing the model's predictive capabilities across various
conditions, like different container geometries or different fluid properties
such as viscosity. All learned from data. But... where is the training data?

So, over the next sections we'll guide you through our recipe of generating synthetic data at scale (and affordably) using our Inductiva API to closely mimic the training data used by the authors of the paper. We will start by setting up a simulation “base case” using SPlisHSPlasH, a popular open-source SPH simulator used by the authors and that is already integrated in our API.
