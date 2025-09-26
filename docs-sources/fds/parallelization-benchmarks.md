## Parallelization Benchmarks
Fire Dynamics Simulator (FDS) simulations benefit greatly from parallelization. FDS supports both MPI and OpenMP parallelization methods:
- **MPI (Message Passing Interface)**: Used to divide the simulation domain into multiple meshes, with each mesh processed in parallel by a separate MPI process.
- **OpenMP**: Enables shared-memory parallelism within each mesh, using multiple threads.

This benchmarks page explores the **performance scaling of FDS simulations** using both MPI and OpenMP. All simulations were run via the **Inductiva API** on **Google Cloud Platform (GCP)** using `c4-standard` machines with **hyperthreading enabled**.

> The "Time to Beat" column in the MPI benchmark table shows reference runtimes obtained from the [official FDS repository](https://github.com/firemodels/out/tree/master/MPI_Scaling_Tests), serving as a baseline to compare against the Inductiva cloud-based simulations.

## MPI Benchmark
To demonstrate the impact of MPI-based parallelization, we replicated the [MPI Strong Scaling benchmark](https://github.com/firemodels/fds/tree/master/Validation/MPI_Scaling_Tests), designed to measure how effectively simulation time decreases as more MPI processes are used.

The folder `FDS_Input_Files` contains simple input cases running for 100 time step. Each case uses a different number of meshes:
- **N=001** → single mesh, run with 1 MPI process (file `strong_scaling_test_001.fds`)
- **N=008** → 8 meshes, run with 8 MPI processes (file `strong_scaling_test_008.fds`)
- **N=016** → 16 meshes, run with 16 MPI processes (file `strong_scaling_test_016.fds`)

The total number of grid cells is kept constant across all cases. Ideally, increasing the number of MPI processes (and hence meshes) should reduce the simulation runtime.

Each simulation was run three times, averaging both runtime and cost. The c4-standard machine with least vCPUs capable of fitting each simulation was selected.

Below are the results for each problem size. The rightmost column shows the corresponding machine cost.

|    MPI Processes |    MPI Slots | Machine Type    |    Avg Time (s) | Avg Cost ($) |    Time to beat (s) |
| ---------------: | -----------: | :-------------- | --------------: | -----------: | ------------------: |
|                1 |            2 | c4-standard-2   |       1360.49   |        0.044 |             1399.00 |
|                8 |            8 | c4-standard-8   |          332.64 |        0.043 |              192.10 |
|               32 |           32 | c4-standard-32  |          116.80 |        0.063 |               62.64 |
|               64 |           96 | c4-standard-96  |           67.04 |        0.117 |               41.54 |
|               96 |           96 | c4-standard-96  |           57.75 |        0.104 |               24.63 |
|              192 |          192 | c4-standard-192 |           37.41 |        0.160 |               14.42 |
|              288 |          288 | c4-standard-288 |           26.39 |        0.167 |                9.80 |

As expected, simulation time decreases as the number of MPI processes increases, demonstrating **effective scaling performance**. Our results compare favorably against the FDS baseline in most configurations.

## OpenMP Benchmark
To demonstrate the effect of OpenMP parallelization, we ran the 8-mesh MPI case with an increasing number of OpenMP threads, while keeping the number of MPI processes fixed at 8. Each case was run on an appropriately sized `c4-standard` machine.

| MPI Processes |    N OMP Threads | Machine Type    |   Avg Time (s) | Avg Cost ($) |
| ------------: | ---------------: | :-------------- | -------------: | -----------: |
|             8 |                1 | c4-standard-8   |         322.60 |        0.042 |
|             8 |                2 | c4-standard-16  |         192.86 |        0.051 |
|             8 |                4 | c4-standard-32  |         139.78 |        0.074 |
|             8 |                6 | c4-standard-48  |         123.16 |        0.099 |
|             8 |               12 | c4-standard-96  |          95.44 |        0.156 |
|             8 |               24 | c4-standard-192 |          84.54 |        0.293 |

Increasing the number of OpenMP threads results in **reduced simulation time**, showcasing the benefits of combining MPI and OpenMP in hybrid parallelization setups. However, cost tends to increase with larger machine sizes, highlighting a trade-off between time and expense.
