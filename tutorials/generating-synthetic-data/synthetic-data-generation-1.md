---
myst:
  html_meta:
    description: "Learn how to generate large-scale synthetic datasets for Physics-ML models using the Inductiva API, starting with our unique data generation recipe."
    keywords: "Inductiva API, Programming, HPC, Simulation, Tutorial, Synthetic Data Generation, Physics-ML, SPH" 
---

# Introduction

How can we train machine learning models to design airfoils for the next
generation of hypersonic planes, optimize the shape of rotor blades for
helicopters destined to fly on Mars, or even predict the interactions between
non-existing molecules and target proteins in drug design, without enough
observational data?

Synthetic data emerges as a solution that addresses the challenge of data
scarcity for training Physics-ML models. Synthetic data generation allows us to
complement or even entirely replace the gaps left by the lack of observational
data with synthetic equivalents. By simulating target systems across a spectrum
of conditions and variations, we can generate diverse and large enough
quantities of synthetic data to train robust Physics-ML models of high accuracy
and strong generalization capabilities. However, we recognize the practical
difficulties and potential expenses associated with generating large and diverse
datasets through simulation, and that's why our Inductiva API stands out as your
go-to tool to make this process more accessible and cost-friendly.

In this tutorial series we'll show you how to use the Inductiva API to generate
synthetic data crucial for training Physics-ML models, at scale. In this
introduction, we're going to let you in on our very own recipe for generating
synthetic data at scale, giving you an overview of the whole process building on
an example from a published study. Over the next steps of this tutorial series,
we'll break it all down to show you how it's done. This guide is a resource for
both machine learning engineers and enthusiasts alike with a thorough
understanding of Physics-ML and simulation software.

## Inductiva's Recipe to Generate Synthetic Data for Physics-ML

Generating synthetic data with the Inductiva API is pretty straightforward and
involves just four steps. From developing an initial simulation model for the
"base case" you wish to study, to running its corresponding simulator across a
number of variations, by the end you'll have a comprehensive and diverse dataset
ready for training ML models.

Our recipe goes as follows:

1. **Set Up the Base Case:** We begin with preparing the configuration files for
   a "base case" simulation model of the system under study. _This step demands
a thorough understanding of the system and simulation software, often requiring
collaboration with specialists in the field. Machine learning engineers will
likely need to partner with domain experts, such as computational chemists or
mechanical engineers, who have experience with specific simulation software. Our
API currently supports a growing range of
[simulators](https://docs.inductiva.ai/en/latest/simulators/overview.html)
covering diverse simulation needs._

2. **Generalize the Base Case:** Next, we transform our "base case" model into a
   template that can generate numerous case variations. _Configuration files
typically model a specific instance of the system with hard-coded values to
describe its parameters. To achieve a dataset rich in diversity and volume, we
need to "generalize" these configuration files in such a way that could be used
to describe relevant variations of the base case. Our [templating
mechanism](https://docs.inductiva.ai/en/latest/explore_api/templating.html)
elegantly replaces fixed values in our configuration files with variables. This
means we can easily tweak these variables through Python scripting later on.
It's a game-changer, allowing us to cover our simulation's parameter space
thoroughly for both exhaustive or random variation runs._

3. **Benchmark Computational Resources:** Next, we benchmark and assess the
   computational performance and costs associated with running variations of our
"base case" across different hardware and at different levels of fidelity.
_Running a large number of high-fidelity simulations can be costly, which is why
we find balancing data resolution against generation time and costs a crucial
step before generating this data at scale. Our API provides all the tools needed
for this benchmarking process with only a few lines of code._

4. **Generate Synthetic Data in Bulk:** Finally, we use our API to fire up
   hundreds of cloud machines, running thousands of variations of our "base
case". Once the simulations complete, we collect the output data for
post-processing. And just like that, we're all set!

Best of all, our API makes it easy to manage our simulations, from listing
resources, to streaming task logs in real-time, even while they're running
remotely, using our [command line interface
(CLI)](https://docs.inductiva.ai/en/latest/cli/cli-overview.html)

We'll get into the details of each of these steps in the coming posts of this
tutorial series. But first, let's take a look at a study done by a group of
researchers who used data generated from Smoothed Particle Hydrodynamics (SPH)
solvers to train a Graph Neural Network (GNN) model.

## Synthetic Data for Learning to Simulate Complex Physics: A Practical Study

Let's take a closer look at a practical application of synthetic datasets
through this comprehensive [study by Sanchez-Gonzalez et
al.](https://arxiv.org/abs/2002.09405). The researchers explore how a Graph
Neural Networks (GNNs) model, trained on synthetic data generated from Smoothed
Particle Hydrodynamics (SPH) solvers, can be applied to simulate fluid dynamics
without learning the laws governing such phenomenon.

The research abstract points out:

> This work uses a state-of-the-art deep learning model to learn the dynamics of
> fluids
and to simulate how they flow without the need for solving the physics equations
that govern the phenomenon. The proposed model is based on a Graph Network
architecture, which takes as input a graph whose nodes represent fluid particles
and edges physical interactions among them. The model performs sequential steps
of computations involving the nodes that are connected by edges to predict how
the fluid flows. (...) The model is trained on data generated with Smoothed
Particle Hydrodynamics, a Lagrangian method in Computational Fluid Dynamics that
simulates fluids as a set of particles.

Smoothed Particle Hydrodynamics (SPH), stands out for its mesh-free approach
that models fluid flows by allowing coordinates to move along with the fluid. It
is widely celebrated for its ability to simulate exceptionally realistic fluid
dynamics, such as those we see in cinematic visual effects.

However, the researchers behind this study consider a much simpler simulation
case: modeling the behavior of a 3D fluid block, made up of many tiny uniformly
dense particles all packed together in a cube shape, suspended in mid-air and
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
dynamics that are not explicitly encoded within the ML model. This challenge
involves generalizing the model's predictive capabilities across various
conditions, like different container geometries or different fluid properties
such as viscosity.

Synthetic data helps us teach machine learning models without needing real-world
data for every scenario, and without directly embedding the physical law for it.
Instead, we use physical simulators that mimic real-world dynamics under various
simulated conditions to produce data that inherently reflects these laws,
helping the model learn on its own. This tutorial series will walk you through
the process of generating such synthetic training datasets to enable your ML
models to learn complex physics with ease.

## Up Next: Setting Up our "Base Case" Simulation

In the [next chapter]({% post_url 2024-03-13-api-synthetic-data-generation-2
%}), we'll guide you through the first step in our recipe of generating
synthetic data at scale and affordably using our Inductiva API. We will teach
you how to create a “base case” using SPlisHSPlasH, a popular open-source SPH
simulator integrated within our API. This is just the beginning of generating
the large, diverse synthetic data you need for your machine learning projects by
making the most of Inductiva's powerful tools.
