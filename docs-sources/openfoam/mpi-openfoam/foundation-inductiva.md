# OpenFOAM-Foundation on Inductiva

To understand the default behavior of OpenFOAM-Foundation when running on Inductiva, let’s start by running the tutorial provided [here](../quick-start.md).

We will run the exact same simulation but with one small change: we will use the default settings for the VM. Your code should look like this:

```python
import inductiva

cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-16",
    spot=True)
```

When you run this simulation, keep in mind that it is originally divided into **6 sub-domains**. Because of this, you will notice that CPU utilization is quite low:

![CPU Usage](../_static/foundation_6_vcpus.png)

This happens because the VM is configured with `threads_by_core=2`, which is the default behavior for virtual machines (see more [here](https://inductiva.ai/guides/how-it-works/machines/hyperthreading)). As a result, your simulation uses only **6 vCPUs out of 16**, which explains the ~37% CPU utilization.

To fully utilize the machine, you have two options:

---

## 1. Utilize All Physical Cores

This approach is straightforward: divide your domain into the same number of **physical cores** and run the simulation.

* The VM provides 16 vCPUs, but only 8 physical cores.
* Running with 8 processes (one per core) results in ~50% CPU utilization:

![CPU Usage](../_static/quick-start/system_metrics_50_2tpc.png)

You can also configure your VM so that the number of available vCPUs matches the number of physical cores (using `threads_by_core=1`). More information is provided [here](https://inductiva.ai/guides/how-it-works/machines/hyperthreading).
This change will make CPU utilization appear as **100%**, due to the fact that the VM will only have 8 vCPUs:

![CPU Usage](../_static/quick-start/system_metrics_100.png)

> **Note**: To clarify. The `c2d-highcpu-16` has 16 vCPUs with the default `threads_per_core=2`. Once we change to `threads_per_core=1` the VM will only have 8 vCPUs, one vCPU per phisical core. Meaning, that using 8 partitions will result in a CPU utilization of 100%.

---

## 2. Utilize All Available vCPUs

Alternatively, you can run the simulation on all 16 vCPUs. To do this, edit the `Allrun` script and replace every instance of `runParallel` with:

```bash
mpirun -np 16 --use-hwthread-cpus <command> -parallel
```

>**Note**:
>
> * You can learn more about MPI on VMs [here](https://inductiva.ai/guides/how-it-works/machines/mpi-on-vms).
> * Don't forget to add the `-parallel` flag after the OpenFOAM command.

This approach results in ~100% CPU usage while utilizing all 16 vCPUs.

*(Insert image showing 100% CPU usage with 16 vCPUs here.)*

---

## Analyzing the Results

There are several possible configurations for running this simulation. Here’s a summary of the execution times and costs for each case:

| Machine Type   | Threads per Core | vCPUs Used | Execution Time | Cost (US$) |
| -------------- | ---------------- | ---------- | -------------- | ---------- |
| c2d-highcpu-16 | 2                | 6          | 2 min 19 sec   | 0.0030     |
| c2d-highcpu-16 | 1                | 8          | 1 min 58 sec   | 0.0026     |
| c2d-highcpu-16 | 2                | 8          | 1 min 56 sec   | 0.0025     |
| c2d-highcpu-16 | 2                | 16         | 1 min 53 sec   | 0.0025     |

From these results:

* Switching from `threads_per_core=2` to `threads_per_core=1` has little impact (differences might be within the margin of normal variation).
* Running this simulation with all 16 vCPUs seems to be the fastest way to do it, but the difference is minimal.

> **Note**:
>
> * This is a very small test case, so the results may not be representative.
> * Each simulation behaves differently; what works best here may not work in other scenarios, and vice versa.

---

Next, we’ll explore the same test on a larger, more realistic simulation.


# Steady-State CFD Simulation of Wind Flow in the Perdigão Region

To better understand atmospheric flow over complex terrain, we conducted a
steady-state CFD simulation of the **Perdigão region in Portugal**. This site
is notable for its two parallel ridges, which generate intricate wind flow
patterns and make it a reference location for atmospheric research.

The simulation was carried out with OpenFOAM’s `simpleFoam` solver on a
structured, terrain-following graded mesh containing **14 million cells**. The
computational domain spanned **30 × 30 × 3 km**, with idealized atmospheric
boundary layer (ABL) conditions applied at the inlet. Turbulence closure was
modeled using the **k–ε model**, and a stepped
**first-order to second-order convection scheme** was employed to ensure better
convergence.

## Analyzing the Results


| Machine Type   | Threads per Core | vCPUs Used | Execution Time | Cost (US$) |
| -------------- | ---------------- | ---------- | -------------- | ---------- |
| c4d-highcpu-96 | 2                | 48          | 9 hrs 20 min  | 15.48     |
| c4d-highcpu-96 | 1                | 48          | 9 hrs 23 min  | 15.58      |
| c4d-highcpu-96 | 2                | 96          | -   | -     |
| c4d-highcpu-48 | 2                | 48         | -   | -     |

From these results:

* Switching from `threads_per_core=2` to `threads_per_core=1` has some negative impact on the performance.
