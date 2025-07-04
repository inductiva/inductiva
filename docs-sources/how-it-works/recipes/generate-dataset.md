# Generate a Dataset
Generating synthetic data with the Inductiva API is a structured process that begins with a base simulation model and scales to thousands of variations, producing a diverse and robust dataset.

The typical workflow includes the following steps:

1. **Set Up the Base Case**
Start by preparing the configuration files for a base case simulation model of the system under study. This step often requires domain expertise and a solid understanding of the simulation software being used.

2. **Generalize the Base Case**
Generalize the configuration files to allow variations of the base case. 
Our [Templating System](https://inductiva.ai/guides/documentation/intro/templating) enables dynamic substitution 
of variables at runtime, making it easy to modify simulation (hyper)parameters through Python. This supports 
both exhaustive and randomized exploration of the configuration space.

3. **Generate Synthetic Data at Scale**
Use the API to launch thousands of simulation variants in parallel on cloud infrastructure. Output data from each simulation is automatically collected and made available for post-processing and downstream use.

Whether you're a Machine Learning engineer or a simulation expert, the Inductiva workflow provides a scalable and efficient solution for generating synthetic datasets to support physics-informed ML models.

---

Ready to dive in? Check out these hands-on tutorials to kickstart your journey:

- [Generate an OpenFOAM Dataset](https://inductiva.ai/guides/openfoam/generate-openfoam-dataset/index)  
- [Create Synthetic Data for Physics-Informed ML](https://inductiva.ai/guides/splishsplash/synthetic-data-for-piml/index)
