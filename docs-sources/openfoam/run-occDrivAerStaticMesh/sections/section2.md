# Run a Single Simulation
Let’s begin by running a single OpenFOAM simulation using the `occDrivAerStaticMesh` case.
All the required input files are already prepared, as outlined in the previous section.

## Code Overview
The Python code required to run a OpenFOAM simulation using the Inductiva API follows a consistent structure. We adapted it for this specific use case, as shown below.

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c2d-highcpu-112",
    data_disk_gb=100,
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

When the simulation is complete, we terminate the machine, download the results and print a summary of the simulation 
as shown as follows.

```
```

As you can see in the “In Progress” line, the part of the timeline that represents the actual execution of the simulation, 
the core computation time of this simulation was approximately --- seconds.

Next, we’ll show you how to scale this same simulation across multiple machines and significantly 
reduce simulation time by a factor of ---!

Stay tuned!