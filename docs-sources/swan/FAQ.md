**Find answers to commonly asked questions about SWAN.**

<br>

# FAQ

## What commands are available to run SWAN?

For this simulator, Inductiva provides a `command` argument that specifies the executable to run. The available commands are as follows:
- **`swanrun`**: This is the default command when no specific command is provided. It is typically used to generate debug files for your simulation. However, `swanrun` does not support MPI clusters, meaning it can only be run on a single machine.
- **`swan.exe`**: This command is compatible with MPI clusters, enabling simulations to run across multiple machines. However, it may occasionally encounter issues when generating debug files.

```python
# Uses swanrun
task = swan.run(input_dir=input_dir,
                sim_config_filename="a11refr.swn",
                command="swanrun",
                on=cloud_machine)

# Uses swan.exe. Note: the simulation file must be called INPUT
task = swan.run(input_dir=input_dir,
                command="swan.exe",
                on=cloud_machine)
```

> **Recommendation**: For simulations on a single machine, Inductiva recommends using
`swanrun`. If you need to run simulations on an MPI cluster, use `swan.exe`.
Also, `swanrun` expects the simulation file in the `file.swn` format, while
`swan.exe` expects the simulation file with the name `INPUT`.

<br>
<br>

Still can't find what you're looking for? [Contact Us](mailto:support@inductiva.ai)
