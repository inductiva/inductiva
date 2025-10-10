# Guide to Bringing Your Own Software
Inductiva is a versatile API platform that simplifies the process of running a wide variety of pre-configured simulation software. With its flexible architecture, Inductiva can also be adapted to run **any scientific software**, making it a powerful tool for researchers and developers.

A key feature of Inductiva is its support for **custom Apptainer images** (formerly Singularity). This allows you to package and upload any software your project requires. Once uploaded, these images can be seamlessly deployed on cloud GPUs, enabling high-performance computing at scale with minimal setup. 

To help you get started and make the most of Inductiva's capabilities, we've created a series of tutorials. These guides are designed to walk you through setting up and running different types of scientific software, giving you the tools and knowledge you need to succeed.

> **Note**: This feature is available exclusively with Inductiva's **Enterprise** 
plan and also included in our **Academia** subscription. For full details, 
please visit our [Pricing](https://inductiva.ai/pricing) page.

---

## ▶️ Watch the Episode
In **Episode 5** of our webinar series *Run It with Inductiva.AI*, we show how to run large-scale, high-performance simulations on the cloud, directly from your laptop, using simple Python scripts.

💡 **What you’ll learn**:
- How to use Inductiva’s Python interface to run scientific computing workflows
- Hands-on simulation examples across selected domains in engineering and science
- Pointers to more examples on scaling simulations and solving related tasks, such as generating synthetic data to train Physics-AI models
- Q&A / Next Steps

<div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%;">
  <iframe src="https://www.youtube.com/embed/BvpUY_b6LcU?si=3yOojtkjmwX_kyB1"
          title="YouTube video player"
          style="position: absolute; top: 0; left: 0; width: 100%; height: 100%; border: 0;"
          allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share"
          allowfullscreen
          referrerpolicy="strict-origin-when-cross-origin">
  </iframe>
</div>

```{banner}
:origin: bring_your_own_software
```

```{toctree}
---
caption: " "
maxdepth: 5
hidden: true
---
how-it-works
run-simulation-with-custom-docker-image
run-sfincs-directly-from-deltares-repository
integrate-your-docker-container
perform-ml-inference/index
run-large-integer-factorization
run-R-stochastic-models
python-based-simulations
```
