# Resource Allocation Options

In this guide, we will explain some of the main features of the Inductiva API
when it comes to making informed decisions about resource allocation. You will
learn how to set up and utilize dedicated resources for running your simulations.
Additionally, we will demonstrate the performance of these resources using a SWASH
simulation.

## Dedicated Resources

To ensure you can access the necessary computational power for your simulations
without the constraints of a shared environment, Inductiva provides the option
to create dedicated pools of Virtual Machine (VM) resources, exclusively reserved
for your use.

There are three types of dedicated computational resources you can launch for your
simulations:

- [**Machine Group**](https://docs.inductiva.ai/en/latest/api_reference/computational_resources/machinegroup_class.html): 
  This consists of homogeneous machines designed to operate individually,
  enabling the distribution of multiple simulations across different machines for
  parallel processing.
  
- [**Elastic Machine Group**](http://docs.inductiva.ai/en/latest/api_reference/computational_resources/elasticgroup_class.html): 
  Similar to Machine Group, these also consist of individual machines. The key
  difference here is the elastic scaling feature, which dynamically adjusts the
  number of machines based on simulation demands, ensuring efficient resource utilization.
  
- [**MPI Cluster**](http://docs.inductiva.ai/en/latest/api_reference/computational_resources/mpicluster_class.html): 
  This setup involves a network of machines configured to work in tandem on a
  single simulation task, distributing the workload across multiple CPUs. This is
  particularly useful for complex simulations that exceed the capabilities of a
  single machine.

## Example: SWASH

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
