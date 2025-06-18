# Run the GRIR443 Benchmark
The GRIR443 benchmark is one of the publicly available cases within the [Quantum ESPRESSO benchmarking suite](https://github.com/QEF/benchmarks/tree/master), widely used to assess the performance of various computational infrastructures. 

This benchmark involves a large-scale self-consistent field (SCF) calculation designed to challenge both memory and 
parallel performance. It models a carbon–iridium complex (C200Ir243) with 2,233,063 G-vectors, eight k-points, and 
FFT grid dimension (180, 180, 192).

According to the [Guide to Running Quantum ESPRESSO (2018)](https://portal.supercomputing.wales/wp-content/uploads/2018/06/Lab_Worksheet_QuantumESPRESSO_SCW_SLURM.pdf?utm_source=chatgpt.com) from the Supercomputing Center of Wales, the GRIR443 benchmark is expected to complete in approximately 22 minutes on 160 CPU cores.

More recently, a performance [report from Fugaku](https://www.hpci-office.jp/documents/appli_software/Fugaku_QE_performance.pdf), the Japanese petascale supercomputer ranked fastest in the world as of June 2020, claimed that the same benchmark completed in just 223.6 seconds when run on 768 AMR-based cores.

In this guide, we will run the GRIR443 benchmark using the Inductiva platform.

## Prerequisites
Download the required files [here](https://github.com/QEF/benchmarks/tree/master/GRIR443) and place them in a folder named `qe_benchmark`. 

## Key Considerations
Running this simulation on Inductiva is as straightforward as running any other Quantum ESPRESSO simulation.

The only aspect to consider is ensuring that there is enough RAM to perform the computation. This means allocating a machine 
with sufficient memory.

Out of the many machines Inductiva makes available, we opted to use the **highmem** variant of the machine type, which comes equipped with 8 GB of RAM per vCPU.

For machines with more than 30 virtual CPUs, it is feasible to run the GRIR443 simulation using the **standard** variant (4 GB 
of RAM per vCPU) or the **highcpu** variant (just 2 GB of RAM per vCPU).

## Run the Simulation
Below is the code required to run GRIR443 with the Inductiva API.

Copy and paste it into a file named `run.py`, save it in the `qe_benchmark` folder, and from that folder, execute it by running:

````
python run.py
````

The simulation will then be sent to a `c3d-highmem-30` virtual machine from Google Cloud, equipped with 30 vCPU and 240GB 
of DDR5 RAM. This machine is equivalent to a high-end laptop, especially due to the amount of RAM available.

```python
import inductiva
from inductiva.commands import MPIConfig, Command

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup( \
    provider="GCP",
    machine_type="c3d-highmem-30",
    spot=True)

mpi_config = MPIConfig( \
    version="4.1.6",
    np=30,
    use_hwthread_cpus=True)

# List of commands to run
commands = [
    Command("pw.x -i grir443.in", mpi_config=mpi_config),
]

# Initialize the Simulator
qe = inductiva.simulators.QuantumEspresso(\
    version="7.4.1")

# Run simulation
task = qe.run( \
    input_dir="GRIR443",
    commands=commands,
    on=cloud_machine)

# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()
task.download_outputs()
task.print_summary()
```

The script will take **~2 hours** to run. In the end, you should see something like:

<task print summary>

That’s it! You have just run the GRIR443 on Inductiva! 

Of course, this run wasn’t as fast as Fugaku’s, nor did it outperform the time reported by the Supercomputing Center of Wales 
in 2018. But we also didn’t use the most powerful machines available on Inductiva!

The magic of Inductiva lies in its flexibility — you can easily scale your simulations to much larger machines, or even multi-node MPI clusters, with just a few minor changes to your Python script.

Let’s explore what happens when we use more powerful machines to run the same simulation.

## Scaling Up GRIR443
Scaling up simply involves updating the `machine_type` parameter when allocating the cloud machine and setting the `np` 
parameter in the MPI configuration to match the number of vCPUs on the selected machine.

The table below shows the execution times across increasingly powerful machines. As the number of vCPUs increases, we 
transition from using the **highmem** variant (8 GB RAM per vCPU) to the **standard** (4 GB RAM per vCPU) and 
**highcpu** variants (2 GB RAM per vCPU). This enables stable pricing while reducing computation time.

| Machine Type       | Execution Time        | Estimated Cost (USD) |
|--------------------|-----------------------|---------------------|
| c3d-highmem-30     | 1 hour, 57 minutes    | 0.91                |
| c3d-standard-60    | 55 minutes, 38 seconds| 0.64                |
| c3d-standard-90    | 41 minutes, 19 seconds| 0.71                |
| c3d-highcpu-180    | 25 minutes, 24 seconds| 0.72                |
| c3d-highcpu-360    | 13 minutes, 20 seconds| 0.76                |

Solid results for less than a dollar!

For just over 70 cents, you can outperform what a national supercomputer achieved back in 2018 by using only a single node.

And here’s the exciting part: with only a modest increase in cost, leveraging a multi-node configuration could potentially allow you to even beat Fugaku’s record.




