# Python Client Usage

> **⚠️ Important**: When using BYOC, you are responsible for all costs incurred by VMs running in your GCP account. Always monitor your GCP console and consider setting up billing alerts to track usage and costs.

## Machine Group Creation

Creating a BYOC machine group is nearly identical to creating a regular Inductiva machine group, with just one key difference: the `byoc=True` parameter. This tells Inductiva to launch the machine in your own GCP account instead of Inductiva's infrastructure.

Create a machine group on your own GCP account:

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
```

**Note**: BYOC machine groups work just like regular machine groups - they will automatically start when you run a simulation if they're not already started. This will create a VM in your GCP account. You can optionally call `mg.start()` explicitly if you want to start the machine before running simulations.

### Configuring Auto-Termination

By default, GCP VMs in BYOC machine groups will automatically terminate after **3 minutes of idle time** to help control costs and reduce the risk of accidentally leaving expensive machines running. This is a critical safety feature that helps prevent unexpected charges from forgotten VMs.

**How it differs from regular machine groups**: Regular Inductiva machine groups are managed entirely by Inductiva and automatically terminate when simulations complete. BYOC machine groups run in your account and require this additional safety mechanism to prevent cost overruns.

You can configure this behavior using the `max_idle_time` parameter:

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

> **⚠️ Warning**: Disabling auto-termination (`max_idle_time=None`) removes the safety mechanism that prevents unexpected charges. You are fully responsible for manually terminating machines. Consider setting up GCP billing alerts and monitoring to track costs if you disable this feature.

## Using the Machine Group

The machine group can be used like any other Inductiva machine group:

```python
# Run a simulation on your BYOC machine group
# The machine will start automatically when you run the simulation
task = simulator.run(input_dir=input_dir, on=mg)
task.wait()
task.download_outputs()
mg.terminate()
```

**Note**: Once you call `simulator.run()`, the simulation will continue running on your GCP VM even if you close your local computer. The machine group will start automatically if it's not already running.

**Machine Visibility**: After creating a BYOC machine group, you can view it in the [Inductiva Console](https://console.inductiva.ai/machine-groups) machines page or using the CLI command `inductiva resources list`. Inductiva will display information about your machine including machine type, start time, status, and usage metrics. This information is reported by the task-runner and helps you monitor your BYOC machines alongside your regular Inductiva resources.

### Retrieving an Existing Machine Group

If you need to access a machine group from a different script or session, you can retrieve it by name:

```python
import inductiva

# Get an existing machine group by name
mg = inductiva.resources.machine_groups.get_by_name("byoc-gcp-machine")

# Use it for simulations
task = simulator.run(input_dir=input_dir, on=mg)
```

## Machine Scaling Workaround

Since BYOC machine groups only support one VM per group (unlike regular machine groups where you can specify `num_machines`), you can create multiple machine groups in a loop to scale up:

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
    # Use each machine group for simulations
    task = simulator.run(input_dir=input_dir, on=mg)
```

```{banner_small}
:origin: launch_machines_on_gcp_sec2
```