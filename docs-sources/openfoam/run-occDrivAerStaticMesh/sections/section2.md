# Run a Single Simulation
First, we will run a single OpenFOAM simulation using the `occDrivAerStaticMesh` 
case. This should be straightforward as all the necessary input files are 
already prepared, as described in the previous section.

## Code Overview
The Python code required to run a OpenFOAM simulation using the Inductiva API follows a consistent structure. We adapted it for this specific use case, as shown below.

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-112",
    spot=True)

# Initialize OpenFOAM stack
openfoam = inductiva.simulators.OpenFOAM(
    version="2412",
    distribution="esi"
    )

task = openfoam.run( \
    input_dir="/Path/to/openfoam-occDrivAerStaticMesh",
    shell_script="./Allrun",
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation as shown as follows.

```
```
