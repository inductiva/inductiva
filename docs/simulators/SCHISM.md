# SCHISM

The [SCHISM](http://ccrm.vims.edu/schismweb/) simulator is an open
source simulator for the simulation of 3D baroclinic circulation
across creek-lake-river-estuary-shelf-ocean scales.

## Example

```python
import inductiva


# Set simulation input directory
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "schism-input-example.zip", True)

# Initialize the Simulator
schism = inductiva.simulators.SCHISM()

# Run simulation with config files in the input directory
task = schism.run(input_dir=input_dir,
                  on=computational_resource,
                  n_vcpus=3,
                  num_scribes=2)

helpers.check_task_status(task)
```