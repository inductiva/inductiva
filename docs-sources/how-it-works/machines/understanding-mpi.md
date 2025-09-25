# Understanding Message Passing Interface (MPI) and Your Machine

When using Inductiva, you may come across something called **MPI**.
MPI, or **Message Passing Interface**, is a system supported by some simulators
that allows your simulation to be split into many smaller tasks that run
simultaneously. This parallel execution can significantly reduce the total
computation time of your simulation.

This guide will help you understand how MPI works and how it interacts with the
machines you use in the cloud.

But before diving into MPI, let’s first understand the basics of cloud machine
configurations.

---

## Understanding Your Cloud Machine Configuration

For the purpose of this guide, we’ll focus only on the **CPU side** of machines
and avoid topics related to memory, as they are less relevant here.

### Machine Naming Scheme

Cloud machine types follow a naming pattern. For example:

* `c2d-highcpu-4`
* `c3d-highcpu-16`

Each name has three parts:

1. **Generation** – (`c2d`, `c3d`): represents the hardware generation.
2. **Configuration** – (`highcpu`): indicates the memory/CPU balance.
3. **vCPU count** – (`4`, `16`): tells you how many virtual CPUs (vCPUs) the machine can have.

For this guide, we’ll focus on the **vCPU count**.

---

### What is a vCPU?

Before explaining vCPUs, let’s quickly look at how a CPU works.

A computer’s **CPU** is its “brain,” where computations happen. The CPU is made
up of multiple **cores**, the individual workers that perform tasks in parallel.

Each physical core can usually handle **two threads** (tasks) simultaneously
using a technique called *hyper-threading*. These threads are exposed as
**virtual CPUs (vCPUs)**.

For example, a CPU with 2 physical cores can provide 4 vCPUs.

![Machine Schema](./_static/machine.png)


So, in the machine `c2d-highcpu-4`:

* 2 physical cores × 2 threads per core = 4 vCPUs

And in `c3d-highcpu-16`:

* 8 physical cores × 2 threads per core = 16 vCPUs

### Changing the Default Behaviour

By default, cloud machines expose **all available vCPUs** (i.e., physical cores × threads per core).
For example:

```python
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-4")
```

This machine type has **2 physical cores** with **2 threads per core**, so the system reports **4 vCPUs**.

However, you can change this behaviour by explicitly setting the number of `threads_per_core`. For example:

```python
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-4",
    threads_per_core=1)
```

With `threads_per_core=1`, each physical core exposes only **one thread**. This means the reported number of vCPUs becomes equal to the number of **physical cores**:

* `c2d-highcpu-4`: 2 physical cores × 1 thread per core = **2 vCPUs**
* `c3d-highcpu-16`: 8 physical cores × 1 thread per core = **8 vCPUs**

This setting can be useful if you want to avoid hyper-threading and run simulations strictly on physical cores.

---

## Understanding MPI

This is where MPI comes into play. MPI divides your simulation into smaller
tasks that vCPUs can work on in parallel. To use MPI, you typically run your
simulator with a command specifying the number of workers (vCPUs) you want.

In most cases, Inductiva (or the simulator) handles MPI configuration for you.
However, in some simulators, you may benefit from setting it up yourself for more control.

---

### Example: Running an OpenFOAM Foundation Simulation

Suppose you have an OpenFOAM-Foundation simulation ready to run. As usual,
you use `runParallel` to execute parts of the simulation in parallel.
Internally, OpenFOAM converts this into an MPI command.

Let’s say you’re running on a `c2d-highcpu-4` machine and you divide your domain
into 4 parts. The command:

```bash
runParallel potentialFoam
```

is internally translated into something like this:

```bash
mpirun -np 4 potentialFoam
```

However, your machine has **2 physical cores** (4 vCPUs). By default, MPI only
“sees” the physical cores, not the vCPUs. This means MPI thinks only 2 slots are
available, and you get an error:

```
There are not enough slots available
```

To fix this, you can tell MPI to also use the vCPUs by adding the flag
`--use-hwthread-cpus`:

```bash
mpirun -np 4 --use-hwthread-cpus potentialFoam
```

This allows your simulation to run across all vCPUs.

> **Note**: If you request a machine with `threads_per_core=1`, the system will only expose physical cores. In this case, the MPI flag `--use-hwthread-cpus` has no effect, since there are no hyper-threaded vCPUs available to use.

👉 **Important detail**: Using `--use-hwthread-cpus` does **not** force MPI to
always schedule tasks on hyper-threaded vCPUs. Instead, MPI first binds one task
to each **physical core**. Only if you request more workers than there are
physical cores will MPI start assigning tasks to hyper-threaded vCPUs. This
ensures physical cores are prioritized, and hyper-threading is only used when
necessary. Example: If you run this simulation with
`mpirun -np 2 --use-hwthread-cpus potentialFoam` each one of the two parts on
run on a diferent phisical core.

In short, the error happens because OpenFOAM-Foundation’s `runParallel` does
not automatically include the `--use-hwthread-cpus` flag when constructing MPI
commands.

> **Note**: In contrast, OpenFOAM-ESI allows you to use `runParallel` across all your vCPUs without needing to add extra flags.

---

## Final Thoughts

You now have a better understanding of how cloud machines are structured and how
MPI interacts with them. This knowledge helps you avoid errors and make better
use of your computing resources.

### Key Takeaways

* **vCPUs = virtual CPUs** created by splitting physical CPU cores into multiple threads.
* **Machine naming** tells you the number of vCPUs available.
* **You can control hyper-threading**: setting `threads_per_core=1` makes vCPUs = physical cores (disables hyper-threading).
* **MPI** enables simulations to run in parallel by dividing tasks among workers.
* By default, MPI only counts **physical cores**, not vCPUs.
* When using `--use-hwthread-cpus`, MPI prioritizes physical cores first and only assigns tasks to hyper-threaded vCPUs if more workers are requested.
* If you set `threads_per_core=1`, the `--use-hwthread-cpus` flag has no effect (only physical cores are visible).
* Some simulators (like OpenFOAM Foundation) may require you to manually adjust MPI commands for optimal performance, while others (like OpenFOAM-ESI) handle this automatically.
