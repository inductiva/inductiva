# LES of Hurricane Boundary Layer 🌪

In this tutorial, we demonstrate how to simulate a simplified hurricane boundary
layer (HBL) using Large Eddy Simulation (LES) on the Inductiva platform. This
setup follows the methodology of Bryan et al. (2017) **[1]**.

The hurricane boundary layer is the lowest portion of the atmosphere in a
tropical cyclone where turbulence, surface friction, and radial inflow are
dominant. This LES setup simulates the internal dynamics of the HBL without full
coupling to large-scale weather systems.

We will also demonstrate Inductiva’s ability to efficiently scale this use case,
starting with a cloud machine equivalent to a typical laptop and then scaling up
to more powerful machines.

## Prerequisites

The required files may be downloaded [here](https://storage.googleapis.com/inductiva-api-demo-files/cm1-les-hurr-bound-layer.zip).
Note that manual download is not necessary since the simulation script will
donwload it.

## Run the Simulation

Below is the code required to run the hurricane boundary layer use case with
the Inductiva API.

The simulation will then be sent to a `c2d-highcpu-16` virtual machine from
Google Cloud, equipped with 16 vCPUs and 32 GB of RAM. This machine is
equivalent to a standard working laptop.

```python
"""LES of Hurricane Bound Layer."""
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True,
)

# Set simulation input directory
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "cm1-les-hurr-bound-layer.zip",
    unzip=True,
)

# Initialize the Simulator
cm1 = inductiva.simulators.CM1(version="21.1")

# Run simulation with config files in the input directory
task = cm1.run(
    input_dir=input_dir,
    sim_config_filename="namelist.input",
    on=cloud_machine,
)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()

task.download_outputs()

task.print_summary()
```

Copy and paste it into a file named `run.py` in your working directory and
execute it by running:

````
python run.py
````

When the simulation is complete, we terminate the machine, download the results
and print a summary of the simulation as shown below.

```
Task status: Success

Timeline:
	Waiting for Input         at 09/07, 14:08:46      0.676 s
	In Queue                  at 09/07, 14:08:47      58.479 s
	Preparing to Compute      at 09/07, 14:09:45      4.285 s
	In Progress               at 09/07, 14:09:49      5772.233 s
		└> 5771.837 s      /opt/openmpi/4.1.6/bin/mpirun --use-hwthread-cpus cm1.exe namelist.input
	Finalizing                at 09/07, 15:46:02      6.8 s
	Success                   at 09/07, 15:46:08      

Data:
	Size of zipped output:    617.07 MB
	Size of unzipped output:  1.21 GB
	Number of output files:   776

Estimated computation cost (US$): 0.15 US$
```

As you can see in the "In Progress" line (the part of the timeline that
represents the actual execution of the simulation), the core computation time
of this simulation was approximately 1 hour and 36 minutes (5772 seconds).

## Scaling Up

One of the strengths of running LES simulations on the Inductiva platform is
the ability to scale your workload across a wide range of high-performance
cloud machines. Whether you need faster results or better cost-efficiency,
scaling up is as simple as adjusting the `machine_type` parameter when
allocating your cloud machine.

We tested the same simulation across several cloud machines with increasing
vCPU counts and next-generation architectures. The results show a clear trend:
more powerful machines can dramatically reduce execution time while maintaining
reasonable costs.

| Machine Type     | Execution Time | Speedup   | Estimated Cost (USD) |
|------------------|----------------|-----------|----------------------|
| c2d-highcpu-16   | 1h 36min       | Reference | $0.14                |
| c4-highcpu-16    | 1h 6min        | 1.45×     | $0.33                |
| c2d-highcpu-112  | 17 min         | 5.65×     | $0.17                |
| c3d-highcpu-360  | 9 min          | 10.67×    | $0.51                |

Upgrading from `c2d-highcpu-16` (a typical entry-level machine) to
`c4-highcpu-16` (a newer generation with the same number of vCPUs) improved
runtime by 45%.

Increasing the number of vCPUs to 112 on `c2d-highcpu-112` slashed the
execution time by over 5×, with only a slight increase in cost—offering. It's a
great balance between speed and cost-efficiency.

For maximum performance, `c3d-highcpu-360` brought the simulation time down to
just 9 minutes, achieving a 10× speedup compared to the baseline. While it’s
the most expensive option tested, it delivers exceptional speed for
time-critical workloads.

These results demonstrate that Inductiva not only supports rapid scaling but
also gives users flexibility to optimize for speed, cost, or a balance of both, depending on the demands of the task. 🚀

If you want to benchmark your own workload in a single script please follow
this [tutorial](https://inductiva.ai/guides/scale-up/benchmark/run-benchmarks).

## References

[1] [Simple hurricane boundary layer (HBL)](https://github.com/george-bryan/CM1/tree/333342b50c85577450868280c2d1cbeff90e2f89/run/config_files/les_HurrBoundLayer)
