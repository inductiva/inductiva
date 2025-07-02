# Set up an MPI Cluster
An **MPI Cluster** is a set of homogeneous machines configured to communicate and collaborate in running a single simulation. By leveraging multiple machines working together, users can execute simulations that require a level of parallelization beyond what a single machine can provide.

Setting up an MPI cluster is straightforward using the `MPICluster` class. You simply specify the type and number of machines to compose the cluster, along with the size of the shared disk storage accessible by all machines.

With this configuration, the Inductiva API initializes the cluster, ensuring all machines can communicate effectively and access the shared storage for reading and writing data.

In the example below, a SWASH simulation is run on an MPICluster. 

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

# Initialize the SWASH simulator and run the simulation
# in your just launched dedicated MPICluster
swash = inductiva.simulators.SWASH()

task = swash.run(
    input_dir=input_dir,
    sim_config_filename="input.sws",
    on=mpi_cluster)

# Wait for the task to finish and download the outputs
task.wait()

# Terminate your dedicated MPICluster at then end of the simulation.
mpi_cluster.terminate()
```

For comparison, the same simulation took 9 minutes and 37 seconds to complete on 
a single `c2-standard-30` machine.

Using the MPI cluster with 120 vCPUs (4x30), the simulation completed in 3 
minutes and 25 seconds, achieving a **2.75× speedup** compared to the single 
machine with 30 vCPUs.

While the time reduction isn’t perfectly linear with the number of vCPUs, 
the improvement remains substantial. For longer simulations, leveraging an 
MPI cluster can drastically reduce execution time from **days to just a few hours**.
