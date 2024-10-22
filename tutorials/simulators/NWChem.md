In this guide, we will walk you through setting up and running NWChem 
simulations, a powerful and versatile computational chemistry software 
built into the Inductiva API. 

We will cover:

- Configuring NWChem for molecular and materials simulations.
- Example code to help you get started with simulations.

# NWChem

NWChem is designed to perform large-scale simulations on modern HPC 
architectures, making it suitable for investigating complex chemical 
phenomena. Its modular architecture allows for ongoing enhancements, 
making it a widely trusted tool for high-accuracy computational chemistry 
in both isolated molecules and extended material systems.

Developed by the Environmental Molecular Sciences Laboratory (EMSL) at 
Pacific Northwest National Laboratory, NWChem offers a range of quantum 
mechanical methods such as Density Functional Theory (DFT), Hartree-Fock, 
and post-Hartree-Fock approaches.

## Example Code

In this example, we run a simple quantum mechanical simulation using NWChem:

```{literalinclude} ../../examples/nwchem/nwchem.py
:language: python
```
