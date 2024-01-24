## SWASH simulator

SWASH is a simulator that solves shallow water equations and is used to simulate 
waves and currents in coastal waters and harbors, long waves in coastal regions 
and tidal inlets, and rapidly varied flows around coastal structures. The simulator 
is configured using a single file with the `.sws` extension, and additional files 
containing information about the domain and the ocean floor, such as a bathymetry 
file with a `.bot` extension, are necessary for the simulation to run.

### Example

```python
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-input-example.zip", unzip=True)

# Initialize the Simulator
swash = inductiva.simulators.SWASH()

# Run simulation with config files in the input directory
task = swash.run(input_dir=input_dir, 
                 sim_config_filename="input.sws")

task.wait()
task.download_outputs()
```