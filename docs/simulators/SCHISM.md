# SCHISM

The [SCHISM](https://ccrm.vims.edu/schismweb/) simulator is an open
source simulator for the simulation of 3D baroclinic circulation
across creek-lake-river-estuary-shelf-ocean scales.

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
    "schism-input-example.zip", True)

# Initialize the Simulator
schism = inductiva.simulators.SCHISM()

# Run simulation with config files in the input directory
task = schism.run(input_dir=input_dir,
                  n_vcpus=3,
                  num_scribes=2,
                  on=machine_group)

task.wait()
task.download_outputs()
```

## Benchmarks

The following benchmarks are currently available for SCHISM:

* [Test_Inun_NWaves_2D](https://benchmarks.inductiva.ai/SCHISM/schism/):
The `Test_Inun_NWaves_2D` example from the SCHISM official test suite.
* [Test_Inun_NWaves_3D](https://benchmarks.inductiva.ai/SCHISM/schism_3d/):
The `Test_Inun_NWaves_3D` example from the SCHISM official test suite.

## What to read next

If you are interested in SCHISM, you may also be interested in checking the
following related simulators that are also avaiable via Inductiva API:

* [Reef3D](Reef3D.md)
* [SWAN](SWAN.md)
* [SWASH](SWASH.md)
* [XBeach](XBeach.md)