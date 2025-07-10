# Preparation for MPI Cluster Execution
In the previous section, we ran the simulation using the `Allrun` shell script — the standard approach for executing 
OpenFOAM simulations. This script wraps the entire simulation workflow into a single command.

However, when scaling the simulation across **multiple machines using an MPI cluster**, we need to take a different approach. 
Rather than relying on `Allrun`, we'll define each simulation step as a separate command. This enables configuring the command 
for parallel execution, fully leveraging the power of the MPI cluster. Don’t worry — this process is handled automatically behind 
the scenes, so there’s no need for manual intervention.

## Setting Up an MPI Cluster
The first step is to update the resource allocation from `MachineGroup` to `MPICluster`, as follows:

```diff
-cloud_machine = inductiva.resources.MachineGroup(
+cloud_machine = inductiva.resources.MPICluster(
    provider="GCP",
    machine_type="c3d-highcpu-180",
+    num_machines=2,
    spot=True)
```

## Defining Your Command List
Next, convert the commands inside the `Allrun` script into a list of standalone commands.

⚠️ **Important Notes**: 
* Each command must be self-contained. For example, directory changes (`cd`) and file operations (`cp`) cannot be split into separate commands. Instead, combine them into a single command, such as:

```bash
cp system/controlDict.noWrite system/controlDict
```

* Special shell characters like `<`, `>`, `|`, and `&` are not allowed, so piping output between commands or redirecting to files is not supported. However, all outputs will be automatically saved to `stdout.txt` and `stderr.txt`.

Here is the result of converting the `Allrun` script into a list of commands:

```python
simulation_commands = [
    "cp system/controlDict.noWrite system/controlDict",
    "cp system/fvSolution.fixedIter system/fvSolution",
    "decomposePar -constant",
    "restore0Dir -processor",
    "renumberMesh -constant -overwrite -parallel",
    "potentialFoam -initialiseUBCs -parallel",
    "applyBoundaryLayer -ybl '0.0450244' -parallel",
    "simpleFoam -parallel",
]
```

## Full Script Using an MPI Cluster
Below is the updated Python script for running the simulation across two machines using an MPI cluster.

```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MPICluster( \
    provider="GCP",
    machine_type="c3d-highcpu-180",
    num_machines=2,
    data_disk_gb=100,
    spot=True)

# Initialize OpenFOAM stack
openfoam = inductiva.simulators.OpenFOAM(
    version="2412",
    distribution="esi")

simulation_commands = [
    "cp system/controlDict.noWrite system/controlDict",
    "cp system/fvSolution.fixedIter system/fvSolution",
    "decomposePar -constant",
    "restore0Dir -processor",
    "renumberMesh -constant -overwrite -parallel",
    "potentialFoam -initialiseUBCs -parallel",
    "applyBoundaryLayer -ybl '0.0450244' -parallel",
    "simpleFoam -parallel",
]

task = openfoam.run( \
    input_dir="/Path/to/openfoam-occDrivAerStaticMesh",
    commands=simulation_commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

Now that we’ve adapted the workflow for an MPI cluster, let’s explore the **performance improvements** we can achieve by scaling 
from one machine to two or even four machines.

