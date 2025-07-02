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

When the simulation is complete, we terminate the machine, download the results
and print a summary of the simulation as shown as follows.

```
```

As you can see in the “In Progress” line, the part of the timeline that
represents the actual execution of the simulation, the core computation time of
this simulation was approximately 3.5 seconds.

Stay with us to see how you can scale this same simulation into multiple machines
and reducing the simulation time by X times! 
