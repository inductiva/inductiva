# The Inductiva Guide to Running Your Own Software
Inductiva is a versatile API platform that simplifies the process of running a wide variety of pre-configured simulation software. With its flexible architecture, Inductiva can also be adapted to run **any scientific software**, making it a powerful tool for researchers and developers.

A key feature of Inductiva is its support for **custom Apptainer images** (formerly Singularity). This allows you to package and upload any software your project requires. Once uploaded, these images can be seamlessly deployed on cloud GPUs,enabling high-performance computing at scale with minimal setup. 

To help you get started and make the most of Inductiva's capabilities, we've created a series of tutorials. These guides are designed to walk you through setting up and running different types of scientific software, giving you the tools and knowledge you need to succeed.

> For technical details on how Apptainer images work within the Inductiva platform, please refer to our [documentation](https://website-staging.inductiva.ai/guides/documentation/intro/private-dockers).

```{toctree}
---
caption: " "
maxdepth: 3
hidden: true
---
run-simulation-with-custom-docker-image
integrate-your-docker-container
perform-ml-inference/index
```
