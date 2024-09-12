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


```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-4",
    num_machines=1,
    data_disk_gb=10)
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "amr-wind-input-example.zip", unzip=True)

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind()

# Run simulation with config files in the input directory
task = amr_wind.run(input_dir=input_dir, 
                 sim_config_filename="abl_amd_wenoz.inp",
                 on=machine_group,
                 n_vcpus=4)

task.wait()
task.download_outputs()

machine_group.terminate()

```

## What to read next

If you are interested in AMR-Wind, you may also be interested in checking the
following related simulators that are also avaiable via Inductiva API:

* [OpenFAST](OpenFAST.md)
* [OpenFOAM](OpenFOAM.md)
