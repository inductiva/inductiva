# Scaling Even Further
In the previous sections of this tutorial, we scaled our simulation from a single machine to a basic MPI cluster with two nodes. 
Now, we’ll take it a step further by running the simulation on four machines, allowing us to explore how performance improves with 
increased computing resources.

In this section, we’ll compare two scenarios:
* With Hyperthreading **Enabled**: Utilizing all available virtual CPUs (vCPUs)
* With Hyperthreading **Disabled**: Using only the physical cores

## Scaling to Four Machines (With Hyperthreading)
The code to run the simulation on four machines is identical to the one used in the two-machine setup — no changes are needed in the Python script itself.

To scale the simulation to four machines, make the following adjustments before executing the script:

1. **Update** the `MPICluster` configuration:
Set `num_machines=4`.

2. **Modify** the `caseDefinition` file to use **448 vCPUs**, as shown below:

```diff
-nCores              224;              // Number of cores used for simulation
+nCores              448;              // Number of cores used for simulation
decompositionMethod hierarchical;      // Decomposition method
-nHierarchical       (14 4 4);           // Coefficients for hierarchical decomposition
+nHierarchical       (14 8 4);          // Coefficients for hierarchical decomposition
```

This setup uses all virtual CPUs available across the four machines, which **includes hyperthreading by default**.

Once these updates are made, you can run the simulation exactly as before.

## Scaling to Four Machines (No Hyperthreading)
To disable hyperthreading and run the simulation using only **physical cores**, adjust the `MPICluster`configuration by setting the `use_hwthread_cpus` parameter to `False` and the `np` parameter to `224` (the total number of physical cores across the four machines).

```python
cloud_machine = inductiva.resources.MPICluster( \
    provider="GCP",
    machine_type="c2d-highcpu-112",
    data_disk_gb=100,
    num_machines=4,
    np=224, #Number of processes MPI will use
    use_hwthread_cpus=False, # Disable hyperthreading
    spot=True)
```

Make sure to revert the `caseDefinition` file to match the new configuration:

```diff
nCores              224;              // Number of cores used for simulation
decompositionMethod hierarchical;      // Decomposition method
nHierarchical       (14 4 4);           // Coefficients for hierarchical decomposition
```

This setup ensures that only the 224 physical cores across the four machines are used.

## Reviewing the Results
Now that you've run the simulation in both configurations, you can compare performance and determine the impact of hyperthreading and cluster scaling.