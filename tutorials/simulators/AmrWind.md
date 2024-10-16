# Amr-Wind

[AMR-Wind](https://github.com/Exawind/amr-wind) stands as a robustly parallel,
block-structured adaptive-mesh solver tailored for simulating incompressible
flows in wind turbines and wind farms. Derived from incflo, its codebase
emphasizes wind-related computations. The solver harnesses the power of the
AMReX library, which furnishes essential components such as mesh data
structures, adaptivity features, and linear solvers crucial for tackling the
governing equations.

AMR-Wind has been built using OpenMPI 4.1.2. We're currently operating on 
AMR-Wind version 1.4.0 for CPU, with a GPU version slated for release in the
near future.

## Example

```{literalinclude} ../../examples/amr-wind/amr-wind.py
:language: python
```

```

## What to read next

If you are interested in AMR-Wind, you may also be interested in checking the
following related simulators that are also avaiable via Inductiva API:

* [OpenFAST](OpenFAST.md)
* [OpenFOAM](OpenFOAM.md)
