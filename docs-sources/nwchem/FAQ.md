**Find answers to commonly asked questions about NWChem.**

<br>

# FAQ

## Why is my NWChem simulation running slowly?
Several factors could contribute to your NWChem simulation running slower than expected. 
One common reason is that the simulation might be too large for the machine you're using. If this is the case, try running the 
simulation on a more powerful machine or a different machine type. Additionally, optimizing your simulation setup can improve performance.

### Common Pitfall: Not Setting the `n_vcpus` parameter
When n_vcpus is not explicitly set, NWChem defaults to using **all available cores** on the selected machine.
This can cause performance issues because NWChem employs both OpenMPI and OpenMP for parallelization. If all the cores are utilized without 
proper tuning, you may **overload the machine**.

For example, if you're using a `c2-standard-4` machine (which has 4 CPU cores), and you don't specify the `n_vcpus` parameter, NWChem will:
- Launch 4 MPI processes (1 per core)
- Each MPI process will spawn 1 OpenMP thread

This results in **4 processes + 4 thread = 8 computation units**—which is not ideal
for a 4-core machine.  

### Recommended Solution
To avoid overloading your machine, explicitly limit the number of cores, like so:

```python
task = nwchem.run(
    input_dir="/Path/to/SimulationFiles",
    sim_config_filename="sim.nw",
    n_vcpus=2,
    on=cloud_machine)
```

> **Tip:** A good rule of thumb is to start with `n_vcpus` set to **half** the
number of cores on your machine. From there, you can experiment with the value to optimize performance.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
