## SWAN

SWAN is a simulator for obtaining realistic estimates of wave parameters in coastal
areas, lakes and eastuaries from given wind, sea floor and current conditions.

The simulator is configured using a single file with the `.swn` extension, and
additional files containing information about the domain and the sea floor, and
the conditions shall be saved in an input directory that is passed to the simulator.

### Example

```python
import inductiva

# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swan-input-example.zip", True)

# Initialize the Simulator
swan = inductiva.simulators.SWAN()

# Run simulation with config files in the input directory
task = swash.run(
    input_dir=input_dir, sim_config_filename="a11refr.swn")

# Wait for the simulation to finish and download the results
task.wait()
task.download_outputs()
```