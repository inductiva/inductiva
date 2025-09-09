**Find answers to commonly asked questions about Quantum ESPRESSO.**

<br>

# FAQ

## 1. What commands are used to run Quantum ESPRESSO?
Here is a list of the commands available to run Quantum ESPRESSO, in alphabetical order:

|                      |                     |                    |                   |                      |
|----------------------|---------------------|--------------------|-------------------|----------------------|
| alpha2f              | band_interpolation  | bands              | bse_main          | casino2upf           |
| cell2ibrav           | cp                  | cppp               | d3hess            | dos                  |
| dynmat               | dvscf_q2r           | epsilon            | ev                | fermi_proj           |
| fermi_velocity       | fs                  | fqha               | gww               | gww_fit              |
| head                 | hp                  | ibrav2cell         | initial_state     | kcw                  |
| kcwpp_interp         | kcwpp_sh            | kpoints            | lambda            | ld1                  |
| manycp               | matdyn              | molecularnexafs    | molecularpdos     | neb                  |
| open_grid            | oscdft_et           | oscdft_pp          | path_interpolation| pawplot              |
| ph                   | phcg                | plotband           | plotproj          | plotrho              |
| pp                   | ppacf               | postahc            | postw90           | pw                   |
| pw2bgw               | pw2critic           | pw2gw              | pw2wannier90      | pw4gww               |
| pwcond               | pwi2xsf             | pprism             | plan_avg          | projwfc              |
| q2qstar              | q2r                 | rism1d             | scan_ibrav        | simple               |
| simple_bse           | simple_ip           | spectra_correction | sumpdos           | turbo_davidson       |
| turbo_eels           | turbo_lanczos       | turbo_magnon       | turbo_spectrum    | upfconv              |
| virtual_v2           | wannier90           | wannier_ham        | wannier_plot      | wfck2r               |
| wfdd                 | xspectra            |                     |                   |                      |

<br>


## 2. Is Quantum ESPRESSO using MPI or OpenMP for parallel execution?

Our Quantum ESPRESSO build supports both **MPI** and **OpenMP** parallelization. The parallelization method depends on the binary you choose to run.

* To use **MPI**, run the standard executable, e.g. `pw.x`. Inductiva will take care of setting the MPI configurations automatically for you, based on the machine you selected before.
* To use **OpenMP**, run the executable with the `_openmp` suffix, e.g. `pw_openmp.x`.

This allows you to select the most suitable parallelization strategy for your simulation.

<br>


## 3. How does Inductiva implicitly set MPI configurations for me?

When you run a **Quantum ESPRESSO** simulation using the **Inductiva API**, the MPI settings are **automatically configured** based on the computational resources you choose.
This means you **don’t need to manually set MPI parameters** for each command, Inductiva handles it for you.

For example:

```python
# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True)
```

For this machine, Inductiva will automatically apply the following MPI configurations:

| MPI Parameter           | Value |
| ----------------------- | ----- |
| **np**                  | 16    |
| **use\_hwthread\_cpus** | True  |

* `np` → Number of MPI processes
* `use_hwthread_cpus` → Whether hyperthreading is enabled

These settings will be applied to **all commands that support MPI**, for example:

```python
commands = [
    "pw.x -i Al_local_pseudo.in"
]
```

<br>

## 4. How can I explicitly configure my MPI settings?

When running Quantum ESPRESSO simulations, you can configure **MPI** settings by specifying the `MPIConfig` parameters in your computational resource.

### Example:

```python
# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-4",
    spot=True)

# Configure MPI settings for the machine
cloud_machine.set_mpi_config(
    mpi_version="4.1.6",      # MPI version to use
    np=2,                     # Number of MPI processes
    use_hwthread_cpus=False)   # Whether to use hyperthreading
```

In this example, the specified MPI configuration will automatically be applied to all commands that support MPI.

For instance, a command like:

```python
command = [
    "pw.x -i Al_local_pseudo.in"
]
```

will internally be executed as:

```bash
mpirun -np 2 --use-hwthread-cpus pw.x -i Al_local_pseudo.in
```

This ensures that the MPI settings are consistently applied across your Quantum ESPRESSO runs.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
