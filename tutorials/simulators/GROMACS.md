In this guide, we will cover how to set up and run GROMACS simulations, 
a versatile molecular dynamics simulator built into the Inductiva API.

We will cover:

- Configuring GROMACS simulations with the necessary input files.
- Example code to help you get started with simulations.

# GROMACS

GROMACS is widely used for simulating biochemical molecules like proteins, 
lipids, and nucleic acids, but its speed and efficiency make it suitable 
for a range of non-biological systems, including polymers and fluid dynamics.

A typical GROMACS simulation consists of several stages: preparing the 
molecular structure, energy minimization, running the simulation, and 
analyzing the results. You’ll need multiple input files to configure and 
run your simulation, and specific GROMACS commands will be used to carry 
out each step.

## Example Code

In the following example, we run a basic molecular dynamics simulation 
of water molecules in a small box:

```{literalinclude} ../../examples/gromacs/gromacs.py
:language: python
```
