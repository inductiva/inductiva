# Set up an MPI Cluster
An **MPI Cluster** is a set of homogeneous machines configured to communicate and collaborate in running a single simulation. By leveraging multiple machines working together, users can execute simulations that require a level of parallelization beyond what a single machine can provide.

Setting up an MPI cluster is straightforward using the `MPICluster` class. You simply specify the type and number of machines to compose the cluster, along with the size of the shared disk storage accessible by all machines.

With this configuration, the Inductiva API initializes the cluster, ensuring all machines can communicate effectively and access the shared storage for reading and writing data.

In the example below, an AMR-Wind simulation is run on an MPICluster. For a step-by-step guide, refer to the related [tutorial](https://inductiva.ai/guides/amr-wind/mpi-cluster-tutorial).

```python
"""AMR-Wind example."""
import inductiva

# Allocate a multi-machine MPI cluster on Google Cloud Platform
cloud_machine = inductiva.resources.MPICluster(
    machine_type="c2d-highcpu-112",
    num_machines=4,
    data_disk_gb=200,
    spot=True
)

# Initialize the Simulator
amr_wind = inductiva.simulators.AmrWind(
    version="3.4.1"
)

# Run the simulation
task = amr_wind.run(
    input_dir="/Path/to/SimulationFiles",
    sim_config_filename="abl_neutral.inp",
    on=cloud_machine
)

# Wait for the simulation to finish and download results
task.wait()
cloud_machine.terminate()
task.download_outputs()
```

The same simulation took **2 hours and 24 minutes** to complete on 
a single `c2d-highcpu-112` machine.

Running it on an MPI cluster with 448 vCPUs (4×112) reduced the runtime to **59 minutes**, resulting in a **2.44× speedup** compared to the single machine with 112 vCPUs.

While the time reduction isn’t perfectly linear with the number of vCPUs, 
the improvement remains substantial. For longer simulations, leveraging an 
MPI cluster can drastically reduce execution time from **days to just a few hours**.
