# CaNS

[CaNS](https://github.com/CaNS-World/CaNS) (Canonical Navier-Stokes) is a
simulator for massively-parallel numerical simulations of fluid flows. It aims
at solving any fluid flow of an incompressible, Newtonian fluid that can benefit
from a FFT-based solver for the second-order finite-difference Poisson equation
in a 3D Cartesian grid.

CaNS has been built using OpenMPI 4.1.2 and FFTW 3.3.8. We're currently
operating on CaNS version 2.3.4 for CPU, with a GPU version slated for release
in the near future.

## Example

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup('c2-standard-4')
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cans-input-example.zip", unzip=True)

# Initialize the Simulator
cans = inductiva.simulators.CaNS()

# Run simulation with config files in the input directory
task = cans.run(input_dir=input_dir, 
                sim_config_filename="input.nml",
                on=machine_group,
                n_vcpus=4)

task.wait()
task.download_outputs()

machine_group.terminate()

```

**Closing Notes**: CaNS requires input files to include a designated 'data'
folder for storing simulation outputs.

## What to read next

If you are interested in CaNS, you may also be interested in checking
the following related simulators that are also avaiable via Inductiva API:

* [DualSPHysics](DualSPHysics.md)
* [OpenFOAM](OpenFOAM.md)
* [SPlisHSPlasH](SPlisHSPlasH.md)
