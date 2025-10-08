**Find answers to commonly asked questions about FDS.**

<br>

# FAQ

## 1. Why aren’t live logs showing up in the Web Console?
FDS behaves differently from most simulators: it writes all logs to `stderr`
instead of `stdout`. Since the Web Console only streams `stdout`, no logs will
appear there during the simulation.

To view the logs while the simulation is running, you can:

* Open the `stderr` file from the live storage.
* Use the CLI command `inductiva logs <task_id>` to stream logs in your terminal.

Once the simulation finishes, both `stdout` and `stderr` files will be available
for viewing in the web console.

<br>

## 2. Why can’t I use more `n_vcpus` for my simulation?
FDS parallelizes simulations by assigning each mesh to a separate vCPU. If your simulation defines only a single mesh, it can only use one vCPU. Attempting to use more will result in an error.

For example, the following setup allows only **1 vCPU**:

```fortran
&MESH IJK=10,10,10 XB=-0.3,0.7,-0.4,0.6,0.0,1.0, MULT_ID='mesh' /
```

To utilize **multiple vCPUs**, define multiple meshes, like so:

```fortran
&MESH IJK= 32,32,144, XB=14.0,18.0, 5.25, 9.25, 0.0,19.0 /
&MESH IJK= 72,21, 72, XB= 0.0,18.0, 0.00, 5.25, 0.0,19.0 /
&MESH IJK= 52,18, 72, XB=14.0,27.0, 9.25,13.75, 0.0,19.0 /
&MESH IJK= 56,34, 72, XB= 0.0,14.0, 5.25,13.75, 0.0,19.0 /
&MESH IJK= 36,37, 72, XB=18.0,27.0, 0.00, 9.25, 0.0,19.0 /
```

This setup allows the simulation to run using 5 vCPUs.

<br>

## 3. How does Fire Dynamics Simulator (FDS) use parallelization?
FDS benefits significantly from parallelization and supports two main methods:
- **MPI (Message Passing Interface)**: This method divides the simulation domain into multiple meshes, with each mesh handled by a separate MPI process running in parallel.
- **OpenMP**: This enables shared-memory parallelism within each mesh by using multiple threads, improving performance on multi-core systems.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
