# Preparation for MPI Cluster Execution
In the previous section, we ran the simulation using the `Allrun` shell script, the standard approach for executing 
OpenFOAM simulations. This script wraps the entire simulation workflow into a single command.

When scaling the simulation across **multiple machines in an MPI cluster**, we
take a different approach. Instead of using `Allrun`, each simulation step is
defined as a separate command. This allows us to configure them for parallel
execution and fully utilize the cluster’s capabilities.

The good news? This entire process is automated, no manual setup required.

## Setting Up an MPI Cluster
The first step is to update the resource allocation from `MachineGroup` to an `MPICluster`, as follows:

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


### How Parallel Commands Are Handled
As mentioned previously, the Inductiva API automatically detects parallel commands
by checking for the `-parallel` flag in your command list. When present, it
configures the command to run across the entire MPI cluster without requiring
any manual setup.

Now that you’ve seen how to move from using the `Allrun` script to a sequence of
individual commands, you're ready to take full advantage of an MPI cluster.
Let’s explore the **performance gains** you can achieve by scaling your simulation
from a single machine to multiple machines.


