# Synthetic Data for Physics-Informed Machine Learning
In this tutorial, we’ll show you how to use the **Inductiva API** to generate synthetic data at scale, efficiently and cost-effectively, for training Physics-informed Machine Learning (PIML) models.

We’ll walk you through our own process for large-scale synthetic data generation, using a real-world example from a published study. Whether you're a machine learning engineer or an enthusiastic practitioner with a solid understanding of simulation software and PIML, this guide is designed for you.

Generating synthetic data with the Inductiva API is straightforward and involves just a few key steps — from building a base simulation model to scaling it across thousands of variations to create a diverse and robust dataset.

Here’s our step-by-step approach:

1. **Set Up the Base Case**
We begin by preparing the configuration files for a base case simulation model of the system under study. This step requires a deep understanding of both the system and the simulation software, often involving collaboration with domain experts. ML engineers may need to work closely with specialists familiar with specific simulators. The Inductiva API currently supports a growing range of simulation engines to cover various physical domains.

2. **Generalize the Base Case**
Next, we generalize the configuration files to represent a wide range of variations of the base case. With our [templating system](https://website-staging.inductiva.ai/guides/documentation/intro/templating), variables within the configuration files are dynamically substituted at runtime. This allows for seamless adjustment of simulation (hyper)parameters through Python scripting, enabling efficient exploration of the configuration space using either exhaustive or randomized strategies.

3. **Benchmark Computational Resources**
Running high-fidelity simulations at scale can be computationally expensive. To ensure efficiency, we first benchmark the trade-offs between data resolution, runtime, and cost. The API includes tools to evaluate different configurations by running the base case across various hardware options and fidelity levels — all with just a few lines of Python code.

4. **Generate Synthetic Data at Scale**
Finally, we use our API to spin up hundreds of cloud machines to run thousands of simulation variations of the base case. Once the simulations complete, we collect the output data for post-processing. And just like that, we’re all set!

To illustrate this process in action, we’ll look at a case study where researchers used data generated from **Smoothed Particle Hydrodynamics (SPH)** simulations to train a **Graph Neural Network (GNN)** model.

```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
sections/section4.md
sections/section5.md
```
