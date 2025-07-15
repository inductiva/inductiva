# Spot Machines

**Spot Machines** (also known as preemptible VMs on GCP) are unused cloud machines available at a significant discount — often 60-90% less than standard on-demand machines. They offer a powerful way to reduce simulation costs. This guide explains how they work and how to use them effectively.

| Feature | On-Demand Machines | Spot Machines |
| :--- | :--- | :--- |
| **Cost** | Standard, fixed price | Heavily discounted (60-90% off) |
| **Reliability** | High (guaranteed availability) | Lower (can be preempted) |
| **Best For** | Time-critical tasks, single runs, and workloads that cannot be interrupted. | Batch processing, fault-tolerant jobs, and cost-sensitive, non-urgent tasks. |

## How Inductiva Handles Preemption

The main drawback of Spot Machines is the risk of interruption. However, Inductiva provides an automatic resubmission mechanism for when you use Spot Machines to run your simulations.

When a Spot Machine is preempted, Inductiva's API will:

1. Detect the interruption.
2. Reschedule the simulation.
3. Relaunch the task on a new machine.

This process is seamless and requires no manual intervention, making Spot Machines a practical and highly cost-effective option for many types of simulations.

You can request Spot Machines for any resource type — MachineGroup, ElasticMachineGroup, or MPICluster — by simply setting the `spot=True` argument during initialization. **Enable automatic resubmission** by setting `resubmit_on_preemption=True` in the simulator's `run()` method.

This ensures you get a low-cost machine _and_ that your task will automatically restart on a new machine if the original is preempted.

```python
import inductiva

# Initialize a MachineGroup with one c2-standard-30 Spot Machine
# The `spot=True` flag requests a cheaper, preemptible instance.
machine_group = inductiva.resources.MachineGroup(
    machine_type="c2-standard-30", 
    spot=True
)

machine_group.start()

swash = inductiva.simulators.SWASH()
task = swash.run(on=machine_group,
                resubmit_on_preemption=True)
```

> Note: If `spot=True` is not set, Inductiva will launch a standard on-demand machine.

## When to Use Spot Machines

Use Spot Machines for:
- Large-scale batch processing and parameter sweeps.
- Simulations that are fault-tolerant by design.
- Non-urgent research and development tasks where cost is a primary concern.
- When running a Benchmark

Avoid Spot Machines for:
- Time-critical simulations with tight deadlines.
- Short, single-run tasks where the potential delay from a preemption would be significant.
- Simulations that cannot be easily or cleanly restarted.


| Pros | Cons |
| :--- | :--- |
| **Significant Cost Savings:** Get access to computing power at a fraction of the on-demand price, often with discounts of 60-90%. | **Preemption Risk:** The cloud provider can reclaim the machine at any time, interrupting your simulation. |
| **Ideal for Flexible Workloads:** Excellent for tasks that are not time-critical or that can tolerate interruptions, such as large batch jobs or parameter sweeps. | **Potential Delays:** If a machine is preempted, the task must be restarted, which can lead to longer overall completion times. |
| **Efficient Resource Utilization:** Leverages the cloud provider's unused capacity, making it a highly efficient choice for scalable computing. | **Not for Critical Deadlines:** Unsuitable for urgent simulations with strict deadlines unless managed by a system that handles preemption automatically. |
