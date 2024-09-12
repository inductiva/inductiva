# SWASH

SWASH is a simulator that solves shallow water equations and is used to simulate 
waves and currents in coastal waters and harbors, long waves in coastal regions 
and tidal inlets, and rapidly varied flows around coastal structures. The
simulator  is configured using a single file with the `.sws` extension, and
additional files containing information about the domain and the ocean floor,
such as a bathymetry file with a `.bot` extension, are necessary for the
simulation to run.

## Example

```python
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-input-example.zip", unzip=True)

# Initialize the Simulator
swash = inductiva.simulators.SWASH()
# or alternatively, to use a specific version of SWASH:
# swash = inductiva.simulators.SWASH(version="10.05")
 

# Run simulation with config files in the input directory
task = swash.run(input_dir=input_dir, 
                 sim_config_filename="input.sws")

task.wait()
task.download_outputs()
```

### Versions

We currently support the following versions of SWASH:
- 9.01A (Apr, 2023)
- 10.01A (Apr, 2024)
- 10.05 (May, 2024)

All available versions for this (and other simulators) can be listed
using the `inductiva simulators list` CLI command.

## Inductiva Benchmarks

The following benchmark is currently available for SWASH:

* [Three-Dimensional Currents](https://benchmarks.inductiva.ai/SWASH/SWASH_Currents/):
This benchmark replicates the S1 simulation as outlined in the paper 
"Modeled Three-Dimensional Currents and Eddies on an Alongshore-Variable Barred
Beach", authored by Christine M. Baker, Melissa Moulton, Britt Raubenheimer, 
Steve Elgar, and Nirnimesh Kumar (2021).

## What to read next

If you are interested in SWASH, you may also be interested in checking the
following related simulators that are also avaiable via Inductiva API:

* [Reef3D](Reef3D.md)
* [SCHISM](SCHISM.md)
* [SWAN](SWAN.md)
* [XBeach](XBeach.md)

You may also want to check the following blog post, where we illustrate a 
practical use of SWASH:

 * [Scaling coastal engineering projects with Inductiva API](https://inductiva.ai/blog/article/scaling-coastal-engineering-projects-inductiva-api)
