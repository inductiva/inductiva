# Resources

The Inductiva API provides flexible computational resources for different simulation needs. Whether you're running a single simulation or orchestrating thousands of parallel tasks, our resource management system handles the complexity of infrastructure provisioning, scaling, and monitoring for you.

## Available Resource Types

Inductiva provides three types of computational resources, each optimized for different simulation patterns and workloads.

### Machine Group
A collection of homogeneous machines designed to operate individually,
enabling the distribution of multiple simulations across different machines for
parallel processing.

- **Best for**: Multiple independent simulations, parameter sweeps, Monte Carlo studies
- **Key benefit**: Parallel execution of separate simulation tasks
- **Learn more**: [Machine Group documentation](computational_resources/machinegroup_class.md)

````{eval-rst}
.. seealso::
   The documentation of the `MachineGroup <https://inductiva.ai/guides/api-functions/api/inductiva.resources#inductiva.resources.machine_groups.MachineGroup>`_ class of out Pyhton API
```` 

### Elastic Machine Group
Similar to Machine Group but with dynamic scaling capabilities that automatically adjust the number of machines based on your simulation queue, ensuring efficient resource utilization and cost optimization.

- **Best for**: Variable workloads, batch processing with unpredictable demand
- **Key benefit**: Automatic scaling based on simulation demand
- **Learn more**: [Elastic Machine Group documentation](computational_resources/elasticgroup_class.md)

### MPI Cluster
A network of machines configured to work collaboratively on a single simulation task, distributing the computational workload across multiple CPUs for maximum performance on complex simulations.

- **Best for**: Complex simulations that exceed the capabilities of a single machine
- **Key benefit**: Parallel processing within a single simulation
- **Learn more**: [MPI Cluster documentation](computational_resources/mpicluster_class.md)

## Quick Start Example

To illustrate the performance of running simulations on dedicated resources using
the Inductiva API, we will run [SWASH](https://inductiva.ai/guides/swash)
as an example.

### Running SWASH on Dedicated Resources

Let's run the simulation on dedicated resources specifically set up for this task.

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

Running the simulation on a [dedicated machine group](#dedicated-resources) with
a `c2-standard-30` machine took **9 minutes and 37 seconds** to complete.
