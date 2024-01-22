## SPlisHSPlasH simulator

SPlisHSPlasH is a Smoothed-Particle Hydrodynamics (SPH) simulator that covers a wide range of applications. The simulator is usually configured by a single file with the extension `.json`. This file contains all the information about the simulation, including the geometry, the physical properties of the fluids, the boundary conditions, the numerical parameters, and the output files. Sometimes the configuration can also use extra geometry files.

### Example

```python
import inductiva

input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "splishsplash-input-example.zip", unzip=True)

# Set simulation input directory
splishsplash = inductiva.simulators.SplishSplash()

task = splishsplash.run(input_dir=input_dir,
                     sim_config_filename="config.json")

task.get_output()
```
