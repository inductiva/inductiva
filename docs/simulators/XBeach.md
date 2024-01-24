## XBeach simulator

XBeach is a simulator with a two-dimensional model for wave propagation, sediment 
transport and morphological changes in the nearshore area. The simulator is configured 
with a `params.txt` file that contains grid and bathymetry info, wave input, flow input, 
morphological input, etc. in the form of keyword/value pairs. If a `params.txt` 
cannot be found then XBeach will not run. Other files are used to configure the 
grid and bathymetry profile, like `bed.dep` for example, and other files with extra 
information that can be used inside the `params.txt` to configure the simulator 
further.

We advise to always set the `mpiboundary` argument in the `params.txt` file, 
since we handle automatically the parallelization of the simulation, based on the 
number of cores available in the machine.

### Example

```python
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "xbeach-input-example.zip", unzip=True)

# Initialize the Simulator
xbeach = inductiva.simulators.XBeach()

# Run simulation with config files in the input directory
task = xbeach.run(input_dir=input_dir,
                  sim_config_filename="params.txt")

task.wait()
task.download_outputs()
```