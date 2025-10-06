# Launching Machines on Your Own GCP Account

Inductiva's Bring Your Own Cloud (BYOC) feature allows you to launch and manage compute machines directly on your own Google Cloud Platform (GCP) account. This gives you full control over your compute resources while leveraging Inductiva's task management and simulation capabilities. Your GCP credentials never leave your local machine, ensuring maximum security and privacy.

## Overview

BYOC enables you to:
- **Reduce costs** by using your own GCP credits and billing account
- **Maintain control** over your compute infrastructure and security policies
- **Scale seamlessly** while keeping compute resources within your organization
- **Comply** with organizational policies that require compute to remain in your cloud account
- **Run simulations seamlessly** with the same experience as Inductiva's managed infrastructure

## What You'll Learn

In this guide, you'll learn how to:
- Set up and configure Google Cloud CLI for BYOC
- Launch machine groups on your own GCP account using the Python client and command line interface
- Configure machine auto-termination
- Use BYOC machine groups for running simulations
- Troubleshoot issues and understand current limitations

## How It Works

When you create a BYOC machine group, Inductiva creates a VM instance in your GCP account that runs a **task-runner container**. This container acts as a bridge between your GCP infrastructure and Inductiva's backend services.

Here's the flow:

1. **VM Creation**: Inductiva creates a VM in your GCP project with the specified machine type and configuration
2. **Task-Runner Deployment**: The VM automatically launches a task-runner container that connects to Inductiva's backend
3. **Seamless Integration**: The task-runner handles communication with Inductiva's API and manages simulation execution
4. **Resource Management**: You maintain control over the VM lifecycle while Inductiva manages the simulation orchestration

**Security Note**: Your GCP credentials never leave your local machine. All VM creation and management operations are performed locally using your authenticated gcloud CLI, ensuring your credentials remain secure and private.

For more detailed information about the task-runner, see the [Task-Runner Guide](../use-local-task-runner/index.md).


## Prerequisites

Before you begin, ensure you have:

- **A Google Cloud Platform (GCP) account** with a valid billing account attached
- **A GCP project** with the Compute Engine API enabled
- **Sufficient permissions** in your GCP project (Compute Instance Admin role recommended)
- **Adequate quotas** for the machine types you plan to use


## Installing and Configuring Google Cloud CLI

### Install Google Cloud CLI

First, install the Google Cloud CLI on your system:

**Linux/macOS:**
```bash
# Download and install
curl https://sdk.cloud.google.com | bash

# Restart your shell or run:
exec -l $SHELL
```

**Windows:**
Download and run the installer from [Google Cloud CLI](https://cloud.google.com/sdk/docs/install)

**Alternative installation methods:**
- See the [official installation guide](https://cloud.google.com/sdk/docs/install) for all options

### Authenticate with Google Cloud

After installation, authenticate with your Google Cloud account:

```bash
# Initialize gcloud
gcloud init
```

The `gcloud init` command will:
1. Open a browser window for authentication
2. Set up your default project
3. Configure your gcloud CLI settings

### Verify Your Setup

Verify that your authentication is working correctly:

```bash
# Check current account
gcloud auth list

# Check current project
gcloud config get-value project
```


### Enable Required APIs

Make sure the Compute Engine API is enabled in your GCP project:

```bash
# Enable Compute Engine API (required for launching VMs)
gcloud services enable compute.googleapis.com
```

You can also enable this through the [GCP Console](https://console.cloud.google.com/apis/library).

## Python Client Usage

### Machine Group Creation

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

### Configuring Auto-Termination

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


### Using the Machine Group

After launching, the machine group can be used like any other Inductiva machine group:

```python
# Run a simulation on your BYOC machine group
task = simulator.run(input_dir=input_dir, on=mg)
task.wait()
task.download_outputs()
mg.terminate()
```

**Note**: Once you call `simulator.run()`, the simulation will continue running on your GCP VM even if you close your local computer.

## Command Line Usage

In addition to the Python client, you can also launch BYOC machine groups directly from the command line using the Inductiva CLI.

### View Available Options

You can view all available flags for the task-runner launch command:

```bash
inductiva task-runner launch -h
```

### Launch a Machine on Your GCP Account

Launch a machine on your own GCP account using the CLI:

```bash
inductiva task-runner launch byoc-gcp-machine \
    --provider gcp \
    --machine-type c2d-highcpu-8 \
    --spot
```

### Access the Machine Group in Python

After launching via CLI, you can access the machine group in Python:

```python
import inductiva

# Get the machine group by name
mg = inductiva.resources.machine_groups.get_by_name('byoc-gcp-machine')

# Use it for simulations
task = simulator.run(input_dir=input_dir, on=mg)
```


## Current Limitations (Alpha Version)

Please be aware of the following limitations in the current alpha release:

### Machine Scaling
- **Single machine per group**: BYOC machine groups only support one VM per group (unlike regular machine groups where you can specify `num_vms`)
- **Workaround**: Create multiple machine groups in a loop to scale up:

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
    ...
    task = simulator.run(input_dir=input_dir, on=mg)
```

### Cost Calculation
- **Cost reporting**: Currently always reports 0 in cost calculations

### Storage
- **Storage support**: Only Inductiva storage is currently supported
- **GCP storage integration**: Coming soon in future releases

## Troubleshooting

### Authentication Errors
**Solution**: Ensure your GCP credentials are properly configured:
- Run `gcloud init` to set up authentication and project configuration
- Verify with `gcloud auth list` that you're authenticated
- Check that your default project is set with `gcloud config get-value project`
- If using service accounts, ensure the service account key is properly configured

### Insufficient Permissions
**Solution**: Ensure your Google account has the necessary IAM role:
- Compute Instance Admin (v1) role
- Check your project's IAM settings in the [GCP Console](https://console.cloud.google.com/iam-admin/iam)

### Insufficient Quotas
**Solution**: 
- Check your GCP project quotas in the [Quotas page](https://console.cloud.google.com/iam-admin/quotas)
- Request quota increases for the specific machine types you need
- Consider using different machine types or regions with available capacity

### Billing Issues
**Solution**: Ensure your GCP project has a valid billing account attached.


```{banner_small}
:origin: launch_machines_on_gcp
```

## Next Steps

- Learn about [exporting files to AWS S3](export-files-to-aws/index.md)
- Explore [bringing your own software](../bring-your-own-software/index.md)
- Check out [local task runner](../use-local-task-runner/index.md) for local development
