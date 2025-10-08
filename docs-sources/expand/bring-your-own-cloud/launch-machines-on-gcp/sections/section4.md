# Command Line Usage

> **⚠️ Important**: When using BYOC, you are responsible for all costs incurred by VMs running in your GCP account. Always monitor your GCP console and consider setting up billing alerts to track usage and costs.

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

> **Note**: The task may be pending for up to two minutes while the machine group starts up.

## Disk Size Configuration

The `data_disk_gb` parameter is very important for BYOC machine groups. Unlike Inductiva's managed infrastructure, BYOC does not support automatic disk resizing. If your simulation runs out of disk space, it will crash and fail.

When launching a BYOC machine group, you can specify the disk size using the `--data-disk-gb` parameter:

```bash
inductiva task-runner launch byoc-gcp-machine \
    --provider gcp \
    --machine-type c2d-highcpu-8 \
    --data-disk-gb 100 \
    --spot
```

Always estimate your simulation's storage requirements and add a safety margin. Monitor your disk usage during simulations to avoid crashes.

## Advanced Configuration

You can customize various aspects of your BYOC machine groups:

```bash
# Launch bigger machine
inductiva task-runner launch byoc-gcp-machine \
    --provider gcp \
    --machine-type c2d-highcpu-32 \
    --spot

# Set custom auto-termination time (in minutes)
inductiva task-runner launch short-lived-machine \
    --provider gcp \
    --machine-type c2d-highcpu-8 \
    --max-idle-time 5 \
    --spot
```

```{banner_small}
:origin: launch_machines_on_gcp_sec3
```
