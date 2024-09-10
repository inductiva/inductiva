# SWAN

SWAN is a simulator for obtaining realistic estimates of wave parameters in coastal
areas, lakes and eastuaries from given wind, sea floor and current conditions.

The simulator is configured using a single file with the `.swn` extension, and
additional files containing information about the domain and the sea floor, and
the conditions shall be saved in an input directory that is passed to the simulator.

## Example

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-4", num_machines=1, data_disk_gb=10)
machine_group.start()

# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swan-input-example.zip", True)

# Initialize the Simulator
swan = inductiva.simulators.SWAN()

# Run simulation with config files in the input directory
task = swan.run(
    input_dir=input_dir, sim_config_filename="a11refr.swn", on=machine_group)

# Wait for the simulation to finish and download the results
task.wait()
task.download_outputs()
```

Check the [official documentation](https://swanmodel.sourceforge.io/) of SWAN to know 
more about the configuration details specific of the simulator.

## Inductiva Benchmarks

The following benchmark is currently available for SWAN:

* [Ring](https://benchmarks.inductiva.ai/SWAN/ring/): The "Ring" example from 
SWAN's site.

## What to read next

If you are interested in SWAN, you may also be interested in checking the
following related simulators that are also avaiable via Inductiva API:

* [Reef3D](Reef3D.md)
* [SCHISM](SCHISM.md)
* [SWASH](SWASH.md)
* [XBeach](XBeach.md)
