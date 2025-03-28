In this guide, we will walk you through setting up and running Quantum ESPRESSO, 
available as one of the built-in tools via the Inductiva API.

We will cover:

- Setting up Quantum ESPRESSO for use with our API.
- Example code to help you get started with simulations.
- The list of allowed commands.

# Quantum ESPRESSO

[Quantum ESPRESSO](https://www.quantum-espresso.org/) is an open-source software
suite widely used for electronic structure calculations and materials modeling
at the nanoscale. It is based on density functional theory (DFT) and uses
plane-wave basis sets to solve quantum mechanical equations for many-body
systems. The package is highly extensible, enabling simulations of a variety of
material properties, including electronic, vibrational, and magnetic
characteristics. Researchers value it for its flexibility, scalability on
high-performance computing platforms, and its role in advancing quantum
simulations and computational materials science.


We have compiled two versions of Quantum ESPRESSO: one for MPI and the
other for OpenMP. To run the MPI version, simply use the standard command names
(e.g., `pw.x`). For the OpenMP version, append `_openmp` to the command names
(e.g., `pw_openmp.x`). This allows users to choose the most suitable version
based on their needs.

All available commands are listed in the [table below](#list-of-allowed-commands). 

## Example

```{literalinclude} ../../inductiva/tests/test_simulators/quantum_espresso/quantum_espresso.py
:language: python
```

## List of allowed Commands
|                      |                     |                   |                   |                    |
|----------------------|---------------------|-------------------|-------------------|--------------------|
| alpha2f              | dvscf_q2r           | head              | matdyn            | plan_avg           |
| pw                   | rism1d              | turbo_spectrum    | average           | dynmat             |
| hp                   | molecularnexafs     | plotband          | pw2bgw            | scan_ibrav         |
| upfconv              | band_interpolation  | epa               | ibrav2cell        | molecularpdos      |
| plotproj             | pw2critic           | simple            | virtual_v2        | bands              |
| epsilon              | initial_state       | neb               | plotrho           | pw2gt              |
| simple_bse           | wannier90           | bse_main          | ev                | kcw                |
| open_grid            | pmw                 | pw2gw             | simple_ip         | wannier_ham        |
| casino2upf           | fermi_proj          | kcwpp_interp      | oscdft_et         | postahc            |
| pw2wannier90         | spectra_correction  | wannier_plot      | cell2ibrav        | fermi_velocity     |
| kcwpp_sh             | oscdft_pp           | postw90           | pw4gww            | sumpdos            |
| wfck2r               | cp                  | fqha              | kpoints           | path_interpolation |
| pp                   | pwcond              | turbo_davidson    | wfdd              | cppp               |
| fs                   | lambda              | pawplot           | ppacf             | pwi2xsf            |
| turbo_eels           | xspectra            | d3hess            | gww               | ld1                |
| ph                   | pprism              | q2qstar           | turbo_lanczos     | dos                |
| gww_fit              | manycp              | phcg              | projwfc           | q2r                |
| turbo_magnon         |                    |                   |                   |                    |


## What to read next

If you are interested in Quantum ESPRESSO, you may also be interested in checking the
following related simulators that are also avaiable via Inductiva API:

* [GROMACS](GROMACS.md)

