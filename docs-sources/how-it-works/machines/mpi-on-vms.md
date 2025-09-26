# Understanding MPI in VMs

MPI enables running parallel tasks across the vCPUs available in your VM. When executing a simulation or computation with MPI, you typically specify the number of processes (`-np`) corresponding to the number of workers you want.

Example MPI command structure:

```bash
mpirun -np <number_of_processes> <program>
```

### How MPI Sees vCPUs

By default, MPI detects only **physical cores**, not hyper-threaded vCPUs. This can cause issues when requesting more workers than physical cores, because MPI does not automatically account for hyper-threaded threads.

For example, a VM with 4 vCPUs but only 2 physical cores will cause MPI to report `There are not enough slots available` unless explicitly told to use hardware threads.

---

### Using Hyper-Threaded vCPUs with MPI

To enable MPI to use hyper-threaded vCPUs, add the flag `--use-hwthread-cpus`:

```bash
mpirun -np <number_of_processes> --use-hwthread-cpus <program>
```

This tells MPI to include hyper-threaded CPUs as valid execution slots, enabling full use of all vCPUs.

---

#### Important Behavior:

* **Without `--use-hwthread-cpus`**: MPI schedules one task per physical core only.
* **With `--use-hwthread-cpus`**: MPI schedules tasks on all vCPUs, starting with physical cores and then using hyper-threaded cores if needed.

Example:
On a VM with 2 physical cores (4 vCPUs):

* `mpirun -np 2 --use-hwthread-cpus ...` → uses only physical cores.
* `mpirun -np 4 --use-hwthread-cpus ...` → uses all vCPUs.

---

### Disabling Hyper-Threading

If your VM is configured with `threads_per_core=1`, only physical cores are available, and hyper-threading is disabled. In this case:

* The MPI flag `--use-hwthread-cpus` has no effect.
* You can only run as many processes as there are physical cores.

---

## Key Takeaways

* **vCPU** = virtual CPU thread provided by hyper-threading or other virtualization technology.
* **VM naming** gives you vCPU count.
* **Hyper-threading** can be controlled with `threads_per_core`.
* **MPI** uses available cores to run tasks in parallel, but by default, counts only physical cores.
* `--use-hwthread-cpus` allows MPI to use hyper-threaded vCPUs.
* Disabling hyper-threading limits MPI to physical cores.
