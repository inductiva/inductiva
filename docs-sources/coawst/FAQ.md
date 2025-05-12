**Find answers to commonly asked questions about COAWST.**

<br>

# FAQ

## Do I need to compile COAWST for every simulation?

**No**, you don’t need to compile COAWST every time — as long as you're using
the same configuration.
If you're running multiple simulations with the same COAWST configuration, you
can **compile it once**, save the resulting binary, and reuse it in subsequent
runs. This can significantly reduce simulation time and resource usage.

You can do this by following the steps below.

### Step 1: Run COAWST and save the binary

On your first run, run COAWST and save the binary in your outputs:

```python
coawst.run(
    input_dir="/Path/to/SimulationFiles",
    sim_config_filename="sim_file.in",
    build_coawst_script="build_coawst.sh",
    compile_simulator=True,
    cleanup_commands=["cp __COAWST/coawstM ."],
    on=machine_group,
    n_vcpus=96)
```

**Notes**:
- `build_coawst_script`is mandatory if you want to compile COAWST.
- `compile_simulator` is set to True by default, so you can omit it if you want.
- `cleanup_commands` is a list of commands that will be executed after the simulation
  is finished. In this case, we are copying the binary `coawstM` to our simulation
  files. You can add any other command you want to run after the simulation is finished.


### Step 2: Reuse the saved binary in future runs

For subsequent simulations, skip compilation and use the saved binary:

```python
coawst.run(
    input_dir="/Path/to/SimulationFiles",
    sim_config_filename="sim_file.in",
    compile_simulator=False,
    # Specify the binary name in your simulation files
    coawst_bin="coawstM",
    on=machine_group,
    n_vcpus=96)
```
**Notes**:
- `build_coawst_script`is not needed if you are not compiling COAWST.
- `compile_simulator` needs to be set to False, so the simulator will not be compiled again.
- `coawst_bin` is the name of the binary you want to use. In this case, we are using
  the binary we saved in the previous simulation.

This is all you need to do to save time and resources when running multiple simulations
using the same COAWST configuration.

> **Note**: If your binary is not called `coawstM`, update the`cleanup_commands` and `coawst_bin` parameters accordingly.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
