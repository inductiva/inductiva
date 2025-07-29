# Resources

The Inductiva API provides flexible computational resources for different simulation needs. Whether you're running a single simulation or orchestrating thousands of parallel tasks, our resource management system handles the complexity of infrastructure provisioning, scaling, and monitoring for you.

## Available Resource Types

Inductiva provides three types of resources you can launch for your simulations.

### Machine Group
A collection of homogeneous machines designed to operate individually,
enabling the distribution of multiple simulations across different machines for
parallel processing.

- **Best for**: Multiple independent simulations
- **Key benefit**: Parallel execution of separate simulation tasks
- **Learn more**: [Machine Group documentation](computational_resources/machinegroup_class.md)

````{eval-rst}
.. seealso::
   For complete API documentation including all parameters, methods, and configuration options, see the `MachineGroup <https://inductiva.ai/guides/api-functions/api/inductiva.resources#inductiva.resources.machine_groups.MachineGroup>`_ class documentation
```` 

### Elastic Machine Group
Similar to Machine Group but with dynamic scaling capabilities that automatically adjust the number of machines based on your simulation queue, ensuring efficient resource utilization and cost optimization.

- **Best for**: Variable workloads, batch processing with unpredictable demand
- **Key benefit**: Automatic scaling based on simulation demand
- **Learn more**: [Elastic Machine Group documentation](computational_resources/elasticgroup_class.md)

````{eval-rst}
.. seealso::
   For complete API documentation including all parameters, methods, and configuration options, see the `ElasticMachineGroup <https://inductiva.ai/guides/api-functions/api/inductiva.resources#inductiva.resources.machine_groups.ElasticMachineGroup>`_ class documentation
````

### MPI Cluster
A network of machines configured to work collaboratively on a single simulation task, distributing the computational workload across multiple CPUs for maximum performance on complex simulations.

- **Best for**: Complex simulations that exceed the capabilities of a single machine
- **Key benefit**: Parallel processing within a single simulation
- **Learn more**: [MPI Cluster documentation](computational_resources/mpicluster_class.md)

````{eval-rst}
.. seealso::
   For complete API documentation including all parameters, methods, and configuration options, see the `MPICluster <https://inductiva.ai/guides/api-functions/api/inductiva.resources#inductiva.resources.machine_groups.MPICluster>`_ class documentation
````

## Choosing the Right Resource

Choose the resource type that best fits your simulation needs:

| | Machine Group | Elastic Machine Group | MPI Cluster |
|---------|---------------|----------------------|-------------|
| **Scaling** | Fixed number of machines | Dynamic auto-scaling | Fixed cluster size |
| **Simulation Distribution** | One simulation per machine | One simulation per machine | Single simulation across all machines |
| **Best for** | Parallel independent tasks | Variable workloads | Large-scale distributed computing |
| **Cost** | Pay for all machines | Pay for active machines only | Pay for entire cluster |
| **Resource Utilization** | All machines run independently | Number of active machines scale with demand | All machines work together |

## Quick Start Example

To illustrate the performance of running simulations on resources using
the Inductiva API, we will run [SWASH](https://inductiva.ai/guides/swash) Simulator on
`c2-standard-30` Google Cloud Platform (GCP) machines as an example. To keep costs down, we will use a preemptible [Spot Machine](spot-machines.md).

```python
import inductiva

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-resources-example.zip", unzip=True)

# Instantiate a MachineGroup object with 1 preemptible machine of type
# c2-standard-30 and start it immediately
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-30", spot=True)
machine_group.start()

# Initialize the SWASH simulator and run the simulation
# in your just launched dedicated MachineGroup
swash = inductiva.simulators.SWASH()
task = swash.run(input_dir=input_dir,
                 sim_config_filename="input.sws",
                 on=machine_group)

# Wait for the task to finish and download the outputs
task.wait()

# Terminate your dedicated MachineGroup at the end of the simulation.
machine_group.terminate()
```

```{banner_small}
:origin: resources
```

```{toctree}
---
caption: Resources
maxdepth: 2
hidden: true
---
computational_resources/machinegroup_class
computational_resources/elasticgroup_class
computational_resources/mpicluster_class
```