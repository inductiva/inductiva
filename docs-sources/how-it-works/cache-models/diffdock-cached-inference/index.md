# Run Generic Scientific Software on Inductiva
Inductiva is a flexible API platform designed to make it easy to run a wide range of pre-configured 
simulation software. However, with a little preparation, Inductiva can also be used to run any 
scientific software, making it a powerful tool for researchers and developers. The platform supports
the use of **custom Apptainer images** (formerly known as Singularity), which allow users to package and upload any software 
required for their tasks. These images can then be efficiently deployed on cloud GPUs.

To showcase the flexibility of Inductiva, this guide walks you through the process of preparing and running [DiffDock](https://github.com/gcorso/DiffDock) — a well-known and somewhat 
complex **Machine Learning model** for molecular docking.

We chose DiffDock not only because it's widely recognized in the scientific community, but also because it's challenging to set up. By understanding how to leverage Inductiva to run DiffDock, 
you’ll be well-equipped to do the same with the software of your choice. 

## Why Run DiffDock on Inductiva?
Running DiffDock on Inductiva offers significant advantages for large-scale computational tasks:
- **Massive parallelization**: You can run hundreds of DiffDock instances simultaneously on cloud GPUs, allowing you to explore millions of protein-ligand pairs quickly and efficiently.
- **Cost-effective resources**: Because DiffDock requires only inference (rather than model training), it can run on relatively inexpensive GPUs, making it ideal for large-scale screening tasks without high computational costs.

By running DiffDock—or **any software of your choice**—on Inductiva, you can harness the power of cloud computing while ensuring that your workflow scales seamlessly with your needs! ⚡️


```{toctree}
:hidden:
sections/section1.md
sections/section2.md
sections/section3.md
```
