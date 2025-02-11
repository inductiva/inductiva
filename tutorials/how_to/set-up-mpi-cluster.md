# Set up an MPI Cluster

The `MPICluster` is composed by a group of homogeneous machines that are enabled
to communicate with each other and work together to run a single simulation.
With several machines working together, users can launch simulations that require
a higher level of parallelization that a single machine cannot offer.

Setting up an MPI cluster is simple with the
[`MPICluster` class](https://docs.inductiva.ai/en/latest/api_reference/computational_resources/mpicluster_class.html).
Users can configure the type and number of machines they want
to compose the cluster and the size of the disk storage shared between all machines.

That is all that is required to launch the cluster. Inductiva API uses this
configuration to initialize the cluster, making sure that all machines can communicate
with each other and can read/write to the shared disk storage.

In the following example, we run a SWASH simulation in an MPICluster. Recall, that
the same simulation took 9m37s to complete on a
[single machine of `c2-standard-30`](https://tutorials.inductiva.ai/intro_to_api/shared_dedicated_resources.html#swash-on-dedicated-resources):

```python
import inductiva

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
   "https://storage.googleapis.com/inductiva-api-demo-files/"
   "swash-resources-example.zip", unzip=True)

# Instantiate a MPICluster object with 4 machine of type c2-standard-30 and 
# start it immediately. This accounts for 120 vCPUs.
mpi_cluster = inductiva.resources.MPICluster(
   provider="GCP",
   machine_type="c2-standard-30",
   num_machines=4)
mpi_cluster.start()

# Initialize the SWASH simulator and run the simulation
# in your just launched dedicated MPICluster
swash = inductiva.simulators.SWASH()

task = swash.run(input_dir=input_dir,
                sim_config_filename="input.sws",
                on=mpi_cluster)

# Wait for the task to finish and download the outputs
task.wait()

# Terminate your dedicated MPICluster at then end of the simulation.
mpi_cluster.terminate()
```

With the MPI cluster using 120 vCPUs, the simulation took 3m25s to finish, which
is a 2.75 times reduction compared to using a single machine with 30 vCPUs.

Notice that the time reduction is not linear with the number of vCPUs, but it
is still a significant improvement. This speed-up can be significant when running
longer simulations and using a readily available MPI cluster can reduce the
time to obtain results from one day to a few hours.
