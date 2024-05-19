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

input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "splishsplash-input-example.zip", unzip=True)

# Set simulation input directory
splishsplash = inductiva.simulators.SplishSplash()

task = splishsplash.run(input_dir=input_dir,
                        sim_config_filename="config.json")

task.wait()
task.download_outputs()
```

## Inductiva Benchmarks

The following benchmarks are currently available for SPlisHSPlasH:

* [Three-Dimensional Currents](https://benchmarks.inductiva.ai/SWASH/SWASH_Currents/):
This benchmark replicates the S1 simulation as outlined in the paper 
"Modeled Three-Dimensional Currents and Eddies on an Alongshore-Variable Barred
Beach", authored by Christine M. Baker, Melissa Moulton, Britt Raubenheimer, 
Steve Elgar, and Nirnimesh Kumar (2021).

## What to read next

If you are interested in SPlisHSPlasH, you may also be interested in checking
the following related simulators that are also avaiable via Inductiva API:

* [CaNS](CaNS.md)
* [DualSPHysics](DualSPHysics.md)
* [OpenFOAM](OpenFOAM.md)

You may also want read the following tutorial, where we exemplify the creation
of synthetic data for training a Physics-ML model using SPlisHSPlasH:

 * [Generating Synthetic Data for training Physics-ML models](https://tutorials.inductiva.ai/generating-synthetic-data/synthetic-data-generation-1.html)