# Preparing for MPI Cluster Execution

In the previous part of this tutorial, we ran the simulation using the `Allrun`
shell script, which is the standard way to execute OpenFOAM simulations. This
script conveniently automates the entire workflow into a single command.

However, when scaling the simulation across multiple machines using an **MPI cluster**,
we need to adjust our approach. Instead of relying on the `Allrun` script, we
will explicitly define each step as separate commands. This allows us to
configure said command for parallel execution, harnessing the full power of your
MPI cluster. Don’t worry, this process is automated behind the scenes, so you
won’t need to manage it manually.

Keep reading to learn how to adapt your simulation for multi-node execution.

## Switching to an MPI Cluster

The first step is to switch your compute resource from a `MachineGroup` to an
**MPICluster**. This is straightforward, simply update your code like this:

```diff
-cloud_machine = inductiva.resources.MachineGroup(
+cloud_machine = inductiva.resources.MPICluster(
    provider="GCP",
    machine_type="c3d-highcpu-180",
+    num_machines=2,
    spot=True)
```

## Defining Your Command List

Next, convert the commands inside your `Allrun` script into a list of
independent commands. A few important notes:

* Each command must be self-contained. For example, you cannot split directory changes (`cd`) and file operations (`cp`) into separate commands. Instead, combine them into a single command like:

  ```bash
  cp system/controlDict.noWrite system/controlDict
  ```

* Special shell characters like `<`, `>`, `|`, and `&` are not allowed, so you cannot pipe output between commands or redirect to files. Don’t worry, all outputs will be automatically saved into `stdout.txt` and `stderr.txt`.

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

---

Your simulation script should now look like this:

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
    distribution="esi"
    )

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

Keep reading to see the speed ups we can get when moving from a single machine into
2 or even 4 machines.