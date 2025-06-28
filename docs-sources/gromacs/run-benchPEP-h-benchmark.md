# Run the benchPEP-h Benchmark
In this tutorial, we demonstrate how to run the GROMACS benchmark developed by the [Max Planck Institute for Biophysical Chemistry (MPINAT)](https://www.mpinat.mpg.de/en), specifically the [PEP-h variant](https://www.mpinat.mpg.de/grubmueller/bench). 

The PEP-h variant test case simulates a small protein in a water box, reflecting typical usage scenarios. It covers multiple stages of the simulation workflow, such as force calculations, neighbor searching, and particle movement, providing a thorough assessment of computational efficiency. 

This benchmark is both memory and computation-intensive, making it ideal for stressing CPU and memory resources. On the CPUs we are considering, the test typically takes between 5 and 15 minutes to complete. The simulation consists of 12 million atoms, with peptides in water, a 2 fs time step, and constrained hydrogen bonds. It is also compatible with GPUs, making it suitable for evaluating both CPU and GPU based hardware.

In this guide, we will run the benchPEP-h benchmark on the cloud using the Inductiva platform.

## Prerequisites
Download the required file [here](https://www.mpinat.mpg.de/benchPEP-h) and place it in a folder named `inputs`. 

## Key Considerations
Running this bechmark on Inductiva is as straightforward as running any other GROMACS simulation.

The only aspect to consider is ensuring that there is enough RAM to perform the computation. This means allocating a machine with sufficient memory.

We verified that the variant **highcpu** – the one with less memory per vCPU – was not able to handle the benchmark, with the simulations consistently crashing due to lack of RAM. 
Hence, we are left with **standard** and **highmem** machines. However, **highmem** variants – the ones with larger RAM sizes – come with a substantial price increase and no performance gain. 

Therefore, for this benchmark, we opted to use the **standard** variant of the machine types.

## Run the Benchmark on a Single Machine
Below is the code required to run benchPEP-h with the Inductiva API. Copy and paste the code into a file named `run.py` and save it:

```python
import inductiva

commands = [
 "gmx mdrun -s benchPEP-h.tpr -pme cpu -bonded cpu -nb cpu -nsteps 1000"
]

# Allocate cloud machine on Google Cloud Platform
machine_group = inductiva.resources.MachineGroup("n2d-standard-32")

# Initialize the Simulator
gromacs = inductiva.simulators.GROMACS()

# Run simulation
task = gromacs.run(input_dir="input_files/",
 commands=commands,
 on=machine_group)


# Wait for the simulation to finish and download the results
task.wait()
cloud_machine.terminate()
task.download_outputs()
task.print_summary()
```

You should have the following folder structure:

```
- benchPEP-h-benchmark/  
  ├── run.py
  └── inputs
    └── benchPEP-h
```


You can now run it with:
```
python run.py
```

The simulation will then be sent to a `n2d-standard-32` virtual machine from Google Cloud, equipped with 32 vCPU and 128GB of RAM.

The simulation will take around **x minutes** to run. In the end, you should see something like:

```python
TODO
```

That’s it! You have just run benchPEP-h on Inductiva! 

The magic of Inductiva lies in its flexibility — you can easily scale your simulations to much larger machines with just a few minor changes to your Python script.

Let’s explore what happens when we use more powerful machines to run the same simulation.

## Scaling Up benchPEP-h Benchmark
Scaling up simply involves updating the `machine_type` parameter when allocating the cloud machine.

We benchmark on the following cloud machines, all of them are the **standard** variant:

| Machine Type | vCPUs Selected for Benchmarking          |
|--------------|------------------------------------------|
| N2D          | 32, 64, 80, 96, 128, 224                 |
| C2D          | 32, 56, 112                              |
| C3D          | 30, 60, 180                              |

First, we will upload the `benchPEP-h` file to the Inductiva platform so that we can avoid having to upload it with each simulation:

```python
inductiva.storage.upload(
    local_path="inputs/benchPEP-h",
    remote_dir="benchPEP-h")
```

This will upload the `benchPEP-h` file to a folder also named `benchPEP-h` on your inductiva storage.


You can now bechmark on multiple machines with the following code, afterwards we will explain what each section of the code does:

```python
from inductiva import simulators, resources, benchmarks

machine_types = [
    "n2d-standard-32",
    "n2d-standard-64",
    "n2d-standard-80",
    "n2d-standard-96",
    "n2d-standard-128",
    "n2d-standard-224",
    "c2d-standard-32",
    "c2d-standard-56",
    "c2d-standard-112",
    "c3d-standard-30",
    "c3d-standard-60",
    "c3d-standard-180",
]

# Create a benchmark
benchmark = benchmarks.Benchmark(name="benchPEP-h")

# Set default simulation parameters
benchmark.set_default(
    simulator=simulators.GROMACS(),
    input_dir=None,
    remote_aseets=["benchPEP-h"],
    commands = ["gmx mdrun -s benchPEP-h.tpr -pme cpu -bonded cpu -nb cpu -nsteps 1000"]
)

# Add runs on the different machines
for machine_type in machine_types:
    benchmark.add_run(on=resources.MachineGroup(machine_type))

# Run the benchmark
benchmark.run(num_repeats=3)

# Wait for tasks to finish
benchmark.wait()

# Export the benchmark results in CSV format
benchmark.export(fmt="csv", filename="benchmark_results.csv")
```

First, we create a benchmark named `benchPEP-h` and set some defaults for all the tasks of that benchmark, namely the `simulator`, the `input_dir`, the `remote_assets` and the `commands` to run.

You can also see that we no longer upload the `inputs` folder in the `input_dir`. This is because we now use the file we previously uploaded to the inductiva storage with the `remote_assets`.

```python
# Create a benchmark
benchmark = benchmarks.Benchmark(name="benchPEP-h")

# Set default simulation parameters
benchmark.set_default(
    simulator=simulators.GROMACS(),
    input_dir=None,
    remote_aseets=["benchPEP-h"],
    commands = ["gmx mdrun -s benchPEP-h.tpr -pme cpu -bonded cpu -nb cpu -nsteps 1000"]
)
```

Then we add a task for each machine type:
```python
for machine_type in machine_types:
    benchmark.add_run(on=resources.MachineGroup(machine_type))
```

We run each task three times to be able to calculte the average:
```python
benchmark.run(num_repeats=3)
```

Finally, we wait for the benchmark to finish and export the results to csv:

```python
benchmark.wait()

benchmark.export(fmt="csv", filename="benchmark_results.csv")
```