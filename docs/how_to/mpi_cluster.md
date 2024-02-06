# Set up MPI Cluster
In this guide, we extend the computational resources to empower the user ability
to run large scaling simulations on a dedicated MPI Cluster.

## How it works?

TODO: Sergio




### Example

Now that we explained what this class is about, let’s see how to use it to run
a SWASH simulation and understand the gains of using an `MPICluster`.

```python
import inductiva

# Download the input files for the SWASH simulation
input_dir = inductiva.utils.download_from_url(
   "https://storage.googleapis.com/inductiva-api-demo-files/"
   "swash-resources-example.zip", unzip=True)

# Instantiate a MPICluster object with 4 machine of type c2-standard-30 and start it
# immediately. This accounts for 120 vCPUs
mpi_cluster = inductiva.resources.MPICluster(
   machine_type="c2-standard-30", num_machines=4)
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
is still a significant improvement. These speed-up are significant when running
longer simulations and using a readily available MPI cluster can reduce the
time to obtain results from one day to a few hours.


## What to read next
* [Set up an MPI]()
* [Get an overview of the CLI]()