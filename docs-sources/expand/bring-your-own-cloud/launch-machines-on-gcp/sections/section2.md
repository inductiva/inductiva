# Python Client Usage

## Machine Group Creation

Create a machine group on your own GCP account using the `byoc=True` parameter:

```python
import inductiva

# Create a BYOC machine group
mg = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-8",  # GCP machine type to use
    spot=True,                     # Use spot instances for cost savings
    mg_name="byoc-gcp-machine",    # Name for your machine group
    byoc=True,                     # Enable Bring Your Own Cloud
)

# Start the machine group
mg.start()
```

## Configuring Auto-Termination

By default, GCP VMs in BYOC machine groups will automatically terminate after **3 minutes of idle time** to help control costs. You can configure this behavior using the `max_idle_time` parameter:

```python
import inductiva

# Set custom idle time (in minutes)
mg = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-8",
    spot=True,
    mg_name="byoc-gcp-machine",
    byoc=True,
    max_idle_time=10  # Terminate after 10 minutes of idle time
)

# Or disable auto-termination entirely
mg = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-8",
    spot=True,
    mg_name="byoc-gcp-machine",
    byoc=True,
    max_idle_time=None  # No auto-termination
)
```

## Using the Machine Group

After launching, the machine group can be used like any other Inductiva machine group:

```python
# Run a simulation on your BYOC machine group
task = simulator.run(input_dir=input_dir, on=mg)
task.wait()
task.download_outputs()
mg.terminate()
```

**Note**: Once you call `simulator.run()`, the simulation will continue running on your GCP VM even if you close your local computer.

## Machine Scaling Workaround

Since BYOC machine groups only support one VM per group (unlike regular machine groups where you can specify `num_vms`), you can create multiple machine groups in a loop to scale up:

```python
import inductiva

# Create multiple BYOC machine groups for scaling
for i in range(3):
    mg = inductiva.resources.MachineGroup(
        provider="GCP",
        machine_type="c2d-highcpu-8",
        spot=True,
        mg_name=f"byoc-gcp-machine-{i}",
        byoc=True,
    )
    mg.start()
    # Use each machine group for simulations
    task = simulator.run(input_dir=input_dir, on=mg)
```

```{banner_small}
:origin: launch_machines_on_gcp_sec2
```
