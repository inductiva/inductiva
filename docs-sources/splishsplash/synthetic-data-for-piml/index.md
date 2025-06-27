# Synthetic Data for Physics-Informed Machine Learning
Training Physics-Informed Machine Learning (PIML) models often requires large, diverse datasets that reflect real-world physical behavior. However, collecting such data through experiments or high-resolution simulations can be time-consuming and expensive. That’s where **synthetic data generation** becomes essential.

In this tutorial, you'll learn how to generate high-quality synthetic datasets at scale using the **Inductiva API** - a fast and cost-effective solution for supporting the training of PIML models.

We'll walk through a real-world example based on a [published study](https://arxiv.org/abs/2002.09405) by Sánchez-González et al., in which a Graph Neural Network (GNN) is trained on data generated using **Smoothed Particle Hydrodynamics (SPH)** to simulate complex fluid dynamics. The study demonstrates how GNNs can learn to model physical systems directly from data, without relying on traditional numerical solvers.

As summarized in the research abstract:

> This work uses a state-of-the-art deep learning model to learn the dynamics of fluids and to simulate how they flow without the need for solving the physics equations that govern the phenomenon. The proposed model is based on a Graph Network architecture, which takes as input a graph whose nodes represent fluid particles and edges physical interactions among them. The model performs sequential steps of computations involving the nodes that are connected by edges to predict how the fluid flows. (…) The model is trained on data generated with Smoothed Particle Hydrodynamics, a Lagrangian method in Computational Fluid Dynamics that simulates fluids as a set of particles.

So, where does this training data come from?

In the following sections, you’ll learn how to:
* Build a base SPH simulation using **SPlisHSPlasH**, the same open-source solver used in the original study (already integrated with our API).
* Configure simulation parameters to reflect real-world variability.
* Scale up to generate thousands of meaningful dataset variations for model training.

Whether you're a machine learning engineer or a simulation expert, this tutorial offers a practical, scalable workflow for generating your own synthetic datasets to power Physics-informed ML models.

<div style="display: flex; justify-content:center">
<video width=500 loop muted autoplay preload="auto">
<source src="../_static/generating-synthetic-data/viscous_flow.mp4" type="video/mp4">
</video>
</div>

<p align="center"><em>Video 2: Dynamics of a viscous fluid simulated with 26,600 particles performed using SPlisHSPlasH via the Inductiva API</em></p>

```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
```
