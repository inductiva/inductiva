# Command Line Usage

In addition to the Python client, you can also launch BYOC machine groups directly from the command line using the Inductiva CLI.

## View Available Options

You can view all available flags for the task-runner launch command:

```bash
inductiva task-runner launch -h
```

## Launch a Machine on Your GCP Account

Launch a machine on your own GCP account using the CLI:

```bash
inductiva task-runner launch byoc-gcp-machine \
    --provider gcp \
    --machine-type c2d-highcpu-8 \
    --spot
```

## Access the Machine Group in Python

After launching via CLI, you can access the machine group in Python:

```python
import inductiva

# Get the machine group by name
mg = inductiva.resources.machine_groups.get_by_name('byoc-gcp-machine')

# Use it for simulations
task = simulator.run(input_dir=input_dir, on=mg)
```

## Advanced Configuration Options

### Custom Machine Types

You can specify different GCP machine types based on your computational needs:

```bash
# High-memory machine for memory-intensive simulations
inductiva task-runner launch memory-intensive-job \
    --provider gcp \
    --machine-type n1-highmem-8 \
    --spot

# High-CPU machine for CPU-intensive tasks
inductiva task-runner launch cpu-intensive-job \
    --provider gcp \
    --machine-type c2d-highcpu-32 \
    --spot
```

### Regional Configuration

Specify a particular GCP region for your machine:

```bash
inductiva task-runner launch regional-machine \
    --provider gcp \
    --machine-type c2d-highcpu-8 \
    --region us-central1 \
    --spot
```

### Custom Auto-Termination

Set custom idle time for cost control:

```bash
inductiva task-runner launch short-lived-machine \
    --provider gcp \
    --machine-type c2d-highcpu-8 \
    --max-idle-time 5 \
    --spot
```

```{banner_small}
:origin: launch_machines_on_gcp_sec3
```
