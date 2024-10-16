# GROMACS

GROMACS is a versatile simulator to perform molecular dynamics simulations. It 
is primarily designed for biochemical molecules like proteins, lipids and
nucleic acids that have a lot of complicated bonded interactions, but since
GROMACS is extremely fast at calculating the nonbonded interactions (that
usually dominate simulations) many groups are also using it for research on
non-biological systems, e.g. polymers and fluid dynamics.

A single simulation of GROMACS via Inductiva API can comprise several steps - 
e.g., preparing the molecules, minimizing the energy of the system, running the
simulation and post-processing. Hence, to configure a simulation of GROMACS the
user may require several files. Moreover, GROMACS has specific commands to run
certain tasks that already use the files in your input directory. 

## Example

This example runs a simple case of the molecular dynamics of water molecules
inside a small box.

```{literalinclude} ../../examples/gromacs/gromacs.py
:language: python
```
