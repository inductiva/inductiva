
## MPI Strong Scaling Benchmark

Fire Dynamics Simulator (FDS) simulations benefit greatly from parallelization. FDS supports parallelization via MPI and OpenMP. For MPI, the domain must be decomposed in meshes, so that each MPI process processes a different mesh in parallel. Within the same mesh, computations can be parallelized via OpenMP. With Inductiva, you have full control of the number of MPI processes and OpenMP threads used in your FDS simulations.
To showcase parallelization via MPI, we will replicate the MPI Strong Scaling benchmark provided in the [official FDS repository](https://github.com/firemodels/fds/tree/master/Validation/MPI_Scaling_Tests), designed to benchmark the parallelization performance of MPI.
The folder `FDS_Input_Files` contains simple input cases running for 100 time steps, where the domain is divided into different numbers of meshes:
- **N=001** → single mesh, 1 MPI process (file `strong_scaling_test_001.fds`)
- **N=008** → 8 meshes, run with 8 MPI processes (file `strong_scaling_test_008.fds`)
- **N=016** → 16 meshes, run with 16 MPI processes (file `strong_scaling_test_016.fds`)
Because the total number of grid cells remains constant, the runtime should ideally decrease as the number of meshes increases.
The benchmark was run on **c4-standard machines from Google Cloud Platform (GCP)** with **hyperthreading enabled**.
The following is an example script for running a simulation case via the Inductiva API:
```python
import inductiva

# Allocate a machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c4-standard-16",
    spot=True)

# Initialize the Simulator
fds = inductiva.simulators.FDS( \
    version="6.9.1")

# Run simulation
task = fds.run( \
    input_dir="FDS_Input_files",
    sim_config_filename="strong_scaling_test_016.fds",
    n_mpi_processes=16,
    n_omp_threads=1,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

 For each simulation case, the c4-standard machine with least vCPUs capable of fitting the simulation was used.  Each case was run three times, averaging both runtime and cost.
Below are the results obtained for each problem simulation size. The rightmost column represents the cost of the machine

|    MPI Processes |    MPI Slots | Machine Type    |    Avg Time (s) | Avg Cost ($) |    Time to beat (s) |
| ---------------: | -----------: | :-------------- | --------------: | -----------: | ------------------: |
|                1 |            2 | c4-standard-2   |       1360.49   |        0.044 |             1399.00 |
|                8 |            8 | c4-standard-8   |          332.64 |        0.043 |              192.10 |
|               32 |           32 | c4-standard-32  |          116.80 |        0.063 |               62.64 |
|               64 |           96 | c4-standard-96  |           67.04 |        0.117 |               41.54 |
|               96 |           96 | c4-standard-96  |           57.75 |        0.104 |               24.63 |
|              192 |          192 | c4-standard-192 |           37.41 |        0.160 |               14.42 |
|              288 |          288 | c4-standard-288 |           26.39 |        0.167 |                9.80 |

We can observe that there the simulation time decreases as we increase the number of MPI processes.

## Scaling simulations with OpenMP

As mentioned, FDS also supports parallelism via OpenMP. To showcase the speed-ups that can be achieved, we'll run the 8 mesh case from the Strong Scaling benchmark with an increasing number of OpenMP threads, also in `c4-standard` machines.

| MPI Processes |    N OMP Threads | Machine Type    |   Avg Time (s) | Avg Cost ($) |
| ------------: | ---------------: | :-------------- | -------------: | -----------: |
|             8 |                1 | c4-standard-8   |         322.60 |        0.042 |
|             8 |                2 | c4-standard-16  |         192.86 |        0.051 |
|             8 |                4 | c4-standard-32  |         139.78 |        0.074 |
|             8 |                6 | c4-standard-48  |         123.16 |        0.099 |
|             8 |               12 | c4-standard-96  |          95.44 |        0.156 |
|             8 |               24 | c4-standard-192 |          84.54 |        0.293 |
Again, we can see that increasing the number of OpenMP threads results in a lower simulation time.
