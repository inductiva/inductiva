# Installation and Setup

## Required Permissions

To minimize security risks, Inductiva requires only the following minimal GCP permissions:

- **Compute Instance Admin (v1)** role, or a custom role with these specific permissions:
  - `compute.instances.create` - To create VMs in your account
  - `compute.instances.delete` - To delete VMs (including auto-termination)
  - `compute.instances.setMetadata` - To set startup script and configuration

> **Security Note**: Your GCP credentials and API keys are **never sent to Inductiva** and Inductiva **never has access to them**. All GCP operations (VM creation, deletion, metadata management) are performed entirely on your local machine using your authenticated gcloud CLI. Inductiva cannot access your other GCP resources, billing information, or data stored in your account.

## Installing and Configuring Google Cloud CLI

### Installation

Install the Google Cloud CLI on your system:

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

## Testing Your Setup

After completing the GCP configuration, test that everything works together:

```bash
# Test GCP authentication
gcloud auth list

# Test project access
gcloud config get-value project

# Test Compute Engine API access
gcloud compute instances list --limit=1
```

If all commands execute without errors, your setup is ready for BYOC.


```{banner_small}
:origin: launch_machines_on_gcp_sec1b
```
