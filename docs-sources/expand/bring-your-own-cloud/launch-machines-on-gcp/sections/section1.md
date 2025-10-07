# Prerequisites and Setup

## Overview

BYOC enables you to:
- **Reduce costs** by using your own GCP credits and billing account
- **Maintain control** over your compute infrastructure and security policies
- **Scale seamlessly** while keeping compute resources within your organization
- **Comply** with organizational policies that require compute to remain in your cloud account
- **Run simulations seamlessly** with the same experience as Inductiva's managed infrastructure

## How It Works

When you create a BYOC machine group, Inductiva creates a VM instance in your GCP account that runs a **task-runner container**. This container acts as a bridge between your GCP infrastructure and Inductiva's backend services.

Here's the flow:

1. **VM Creation**: Inductiva creates a VM in your GCP project with the specified machine type and configuration
2. **Task-Runner Deployment**: The VM automatically launches a task-runner container that connects to Inductiva's backend
3. **Seamless Integration**: The task-runner handles communication with Inductiva's API and manages simulation execution
4. **Resource Management**: You maintain control over the VM lifecycle while Inductiva manages the simulation orchestration

**Security Note**: Your GCP credentials never leave your local machine. All VM creation and management operations are performed locally using your authenticated gcloud CLI, ensuring your credentials remain secure and private.

For more detailed information about the task-runner, see the [Task-Runner Guide](../../use-local-task-runner/index.md).

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

```{banner_small}
:origin: launch_machines_on_gcp_sec1
```
