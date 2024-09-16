# SPlisHSPlasH

SPlisHSPlasH is a Smoothed-Particle Hydrodynamics (SPH) simulator that covers a 
wide range of applications. The simulator is usually configured by a single file 
with the extension `.json`. This file contains all the information about the 
simulation, including the geometry, the physical properties of the fluids, the 
boundary conditions, the numerical parameters, and the output files. Sometimes
the  configuration can also use extra geometry files.

## Example

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup('c2-standard-4')
machine_group.start()

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "splishsplash-input-example.zip", unzip=True)

# Set simulation input directory
splishsplash = inductiva.simulators.SplishSplash()

task = splishsplash.run(input_dir=input_dir,
                        sim_config_filename="config.json",
                        on=machine_group)

task.wait()
task.download_outputs()

machine_group.terminate()
```

## Inductiva Benchmarks

The following benchmarks are currently available for SPlisHSPlasH:

* [Fluid Cube S](https://benchmarks.inductiva.ai/SPlisHSPlasH/splish_splash/):
This benchmark replicates the example in [this tutorial](https://tutorials.inductiva.ai/generating-synthetic-data/synthetic-data-generation-1.html), using the deafult values.
* [Fluid Cube L](https://benchmarks.inductiva.ai/SPlisHSPlasH/splish_splash/):
This benchmark mirrors the benchmark [Fluid Cube S](https://benchmarks.inductiva.ai/SPlisHSPlasH/splish_splash/) in all aspects except for the particle radius, which has been decreased to 0.0045 and the fluid model, that was doubled in all axis.
* [Fluid Cube M](https://benchmarks.inductiva.ai/SPlisHSPlasH/splish_splash/): This benchmark mirrors the [Fluid Cube S](https://benchmarks.inductiva.ai/SPlisHSPlasH/splish_splash/) benchmark in all aspects except for the particle radius, which has been decreased to 0.0045.

## What to read next

If you are interested in SPlisHSPlasH, you may also be interested in checking
the following related simulators that are also avaiable via Inductiva API:

* [CaNS](CaNS.md)
* [DualSPHysics](DualSPHysics.md)
* [OpenFOAM](OpenFOAM.md)

You may also want read the following tutorial, where we exemplify the creation
of synthetic data for training a Physics-ML model using SPlisHSPlasH:

 * [Generating Synthetic Data for training Physics-ML models](https://tutorials.inductiva.ai/generating-synthetic-data/synthetic-data-generation-1.html)