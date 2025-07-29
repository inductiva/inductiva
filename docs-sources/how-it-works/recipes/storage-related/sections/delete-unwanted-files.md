# Cleaning Up After Simulations

Simulations often generate a large number of temporary or intermediate files
that are only useful during runtime. These files can consume significant
storage, especially when running multiple or large-scale simulations. To help
manage storage usage and reduce associated costs, it is recommended to clean up
unnecessary files once a simulation finishes.

All simulators allow the use of the `on_finish_cleanup` argument to specify
cleanup routines. This argument can be either:
- A path to a shell script
- A list of shell commands

The cleanup routine runs automatically after the simulation ends, making it easy
to keep your workspace tidy and efficient.

## Why It Matters

- **Reduce storage usage**: Prevent accumulation of unused data
- **Minimize costs**: Avoid paying for unnecessary storage when running in the cloud
- **Keep things organized**: Cleaner file structures are easier to debug and manage

## Example: Using a Shell Script

Some simulators, like OpenFOAM, create per-process folders that are only useful
during the run and safe to delete afterward. Here's how to define a cleanup script:

**`clean.sh`**
```bash
#!/bin/bash
echo "Cleaning up process folders..."
rm -rf processor*
```

> **Note**: This clean.sh needs to be passed with your `input_files`.

And then, in your simulation script:

```python
# Run the simulation with a cleanup script
task = openfoam.run( \
    input_dir="/Path/to/your/input_files",
    commands=simulation_commands,
    on_finish_cleanup="clean.sh",
    on=cloud_machine)
```

## Example: Using a List of Shell Commands

If you prefer not to create a separate script, you can pass a list of shell
commands directly:

```python
# Run the simulation with inline cleanup commands
task = openfoam.run( \
    input_dir="/Path/to/your/input_files",
    commands=simulation_commands,
    on_finish_cleanup=[
        "echo 'Cleaning up temporary folders...'",
        "rm -rf postProcessing/temp_data"
    ],
    on=cloud_machine)
```

```{banner_small}
:origin: recipes_delete_unwanted_files
```