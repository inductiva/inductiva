# Run the Simulation on an MPI Cluster
You’re almost ready to scale your simulation across multiple machines, unlocking
the full power of an **MPI Cluster**. 

Running simulations in parallel across a cluster can dramatically reduce computation time by distributing the workload over 
many machines.

## Update Simulation Parameters to Use 360 vCPUs
Currently, the simulation is configured to run on 180 vCPUs. To fully utilize
the MPI Cluster, which provides 360 vCPUs across two machines, you’ll need to
adjust the simulation parameters accordingly.

Open the file `system/include/caseDefinition` and update the number of cores and decomposition settings:

```diff
-nCores              180;              // Number of cores used for simulation
+nCores              360;              // Number of cores used for simulation
decompositionMethod hierarchical;      // Decomposition method
-nHierarchical       (9 5 4);           // Coefficients for hierarchical decomposition
+nHierarchical       (15 6 4);          // Coefficients for hierarchical decomposition
```

These changes ensure that the simulation workload is properly partitioned across all available processors for optimal parallel performance.

## Run your simulation
Below is the updated Python script for running the simulation across two machines using an MPI cluster and the newly
created list of commands.

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

Just like in the single-machine simulation, once the run is complete, the
machines are automatically shut down, the results are downloaded, and a detailed
summary is printed, as shown below.

```
inductiva tasks info 4trnwrg7xscpk6ynzsw7dnxsh
Task status: Success

Timeline:
	Waiting for Input         at 09/07, 10:45:29      0.879 s
	In Queue                  at 09/07, 10:45:30      56.879 s
	Preparing to Compute      at 09/07, 10:46:27      113.184 s
	In Progress               at 09/07, 10:48:20      28246.922 s
		├> 1.08 s          mv caseDefinition system/include/caseDefinition
		├> 1.096 s         cp system/controlDict.noWrite system/controlDict
		├> 1.079 s         cp system/fvSolution.fixedIter system/fvSolution
		├> 129.235 s       decomposePar -constant
		├> 3.09 s          restore0Dir -processor
		├> 63.279 s        /opt/openmpi/4.1.6/bin/mpirun --hostfile /home/task-runner/mpi_hosts --mca btl_tcp_if_include 10.132.0.0/20 --mca oob_tcp_if_include 10.132.0.0/20 --np 360 --use-hwthread-cpus renumberMesh -constant -overwrite -parallel
		├> 45.22 s         /opt/openmpi/4.1.6/bin/mpirun --hostfile /home/task-runner/mpi_hosts --mca btl_tcp_if_include 10.132.0.0/20 --mca oob_tcp_if_include 10.132.0.0/20 --np 360 --use-hwthread-cpus potentialFoam -initialiseUBCs -parallel
		├> 20.108 s        /opt/openmpi/4.1.6/bin/mpirun --hostfile /home/task-runner/mpi_hosts --mca btl_tcp_if_include 10.132.0.0/20 --mca oob_tcp_if_include 10.132.0.0/20 --np 360 --use-hwthread-cpus applyBoundaryLayer -ybl 0.0450244 -parallel
		└> 27974.321 s     /opt/openmpi/4.1.6/bin/mpirun --hostfile /home/task-runner/mpi_hosts --mca btl_tcp_if_include 10.132.0.0/20 --mca oob_tcp_if_include 10.132.0.0/20 --np 360 --use-hwthread-cpus simpleFoam -parallel
	Finalizing                at 09/07, 18:39:07      21.577 s
	Success                   at 09/07, 18:39:28      

Data:
	Size of zipped output:    5.56 GB
	Size of unzipped output:  10.82 GB
	Number of output files:   6517

Estimated computation cost (US$): 26.97 US$

Go to https://console.inductiva.ai/tasks/4trnwrg7xscpk6ynzsw7dnxsh for more details.
```

As shown in the “In Progress” line, all commands, including those running in
parallel, are listed. These commands were executed across **360 vCPUs**,
utilizing the full capacity of our MPI cluster.

With this setup, the total computation time dropped to **7 hours and 53 minutes**,
nearly halving the runtime compared to the single-machine execution.

For this particular simulation, scaling to additional machines offers no further
benefit, as the **communication overhead outweighs** the performance gains from
the extra compute power.

You can learn more about running OpenFOAM on MPI clusters in our [blog post](https://inductiva.ai/blog/article/from-supercomputer-to-cloud-a-new-era-for-openfoam-simulations).


