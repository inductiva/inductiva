# Frequently Asked Questions (FAQ)

## General Questions

### What is BYOC and why should I use it?

BYOC (Bring Your Own Cloud) allows you to launch and manage compute machines directly on your own Google Cloud Platform (GCP) account while using Inductiva's simulation capabilities. This gives you:

- **Cost control**: Use your own GCP credits and billing account
- **Security**: Your GCP credentials never leave your local machine
- **Compliance**: Keep compute resources within your organization
- **Flexibility**: Full control over your compute infrastructure

### Is BYOC secure?

Yes, BYOC is designed with security as a priority:

- **Your GCP credentials never leave your local machine** - they are never sent to Inductiva
- **Inductiva never has access to your GCP credentials** - all GCP operations are performed locally using your authenticated gcloud CLI
- **Minimal permissions required** - only specific compute instance permissions are needed
- **Task-runner isolation** - the container running on your VM acts as a secure bridge

### What does Inductiva see when I use BYOC?

Inductiva can see machine information reported by the task-runner for monitoring and pricing purposes:
- Machine specifications (vCPUs, RAM, disk size)
- Usage metrics (CPU, RAM, disk usage)
- Runtime information (start time, last seen)
- Live logs and file tracking for running simulations

However, Inductiva **cannot** see your GCP credentials, billing information, or access your GCP console.

## Setup and Configuration

### What permissions do I need in GCP?

You need the **Compute Instance Admin (v1)** role, or a custom role with these specific permissions:
- `compute.instances.create` - To create VMs
- `compute.instances.delete` - To delete VMs (including auto-termination)
- `compute.instances.setMetadata` - To set startup script and configuration

### How do I set up the Google Cloud CLI?

1. **Install gcloud CLI**:
   ```bash
   # Linux/macOS
   curl https://sdk.cloud.google.com | bash
   exec -l $SHELL
   ```

2. **Authenticate**:
   ```bash
   gcloud init
   ```

3. **Enable required API**:
   ```bash
   gcloud services enable compute.googleapis.com
   ```

4. **Verify setup**:
   ```bash
   gcloud auth list
   gcloud config get-value project
   ```

### What APIs need to be enabled?

You need to enable the **Compute Engine API** (`compute.googleapis.com`) in your GCP project.

## Machine Management

### How do I create a BYOC machine group?

Simply add `byoc=True` to your machine group creation:

```python
import inductiva

mg = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-8",
    spot=True,
    byoc=True,  # This enables BYOC
)
```

### Can I scale BYOC machine groups?

**Current limitation**: BYOC machine groups only support one VM per group (unlike regular machine groups).

**Workaround**: Create multiple machine groups in a loop:

```python
for i in range(3):
    mg = inductiva.resources.MachineGroup(
        provider="GCP",
        machine_type="c2d-highcpu-8",
        spot=True,
        byoc=True,
    )
    # Use each machine group for simulations
```

### How does auto-termination work?

By default, BYOC VMs automatically terminate after **3 minutes of idle time** to prevent unexpected charges. You can:

- **Customize the timeout**:
  ```python
  mg = inductiva.resources.MachineGroup(
      # ... other parameters ...
      max_idle_time=10  # 10 minutes
  )
  ```

- **Disable auto-termination**:
  ```python
  mg = inductiva.resources.MachineGroup(
      # ... other parameters ...
      max_idle_time=None  # No auto-termination
  )
  ```

⚠️ **Warning**: Disabling auto-termination removes the safety mechanism. You're fully responsible for manually terminating machines.

### How do I configure disk size?

Unlike Inductiva's managed infrastructure, BYOC doesn't support automatic disk resizing. You must specify the disk size upfront:

```python
mg = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-8",
    spot=True,
    byoc=True,
    data_disk_gb=100  # Specify disk size in GB
)
```

**Important**: If your simulation runs out of disk space, it will crash. Always estimate your storage requirements and add a safety margin.

## Cost and Billing

### Who pays for the VMs?

**You are responsible for all costs** incurred by VMs running in your GCP account. This includes:
- VM compute costs
- Storage costs
- Network egress costs
- Any other GCP charges

### How can I monitor costs?

- **Monitor your GCP console** regularly
- **Set up billing alerts** in GCP to track usage
- **Use GCP's cost management tools** to analyze spending
- **Check the Inductiva Console** for machine usage metrics

### Does Inductiva charge for BYOC?

Inductiva doesn't charge for the compute resources (since they run in your account), but standard Inductiva usage fees may still apply for simulation services.

### Why does cost reporting show 0?

This is a current limitation in the alpha version. Cost calculations always report 0 for BYOC machine groups.

## Troubleshooting

### I'm getting authentication errors. What should I do?

1. Run `gcloud init` to set up authentication
2. Verify authentication: `gcloud auth list`
3. Check your default project: `gcloud config get-value project`
4. Ensure you're authenticated with the correct account

### I'm getting permission denied errors. How do I fix this?

1. **Check your GCP IAM permissions** in the [GCP Console](https://console.cloud.google.com/iam-admin/iam)
2. **Ensure you have the required permissions**:
   - `compute.instances.create`
   - `compute.instances.delete`
   - `compute.instances.setMetadata`
3. **Contact your GCP administrator** if you need additional permissions

### My machine isn't starting. What could be wrong?

1. **Check VM status** in the GCP Console
2. **Verify the task-runner container** is running
3. **Check VM logs** in the GCP Console for startup errors
4. **Ensure your project has a valid billing account** attached
5. **Check if you've exceeded quotas** for the machine type

### I'm getting quota exceeded errors. How do I resolve this?

1. **Check your GCP project quotas** in the [Quotas page](https://console.cloud.google.com/iam-admin/quotas)
2. **Request quota increases** for the specific machine types you need
3. **Try different machine types** or regions with available capacity
4. **Contact GCP support** for quota increase requests

### I'm having network connectivity issues. What should I check?

1. **Ensure outbound HTTPS traffic (port 443)** is allowed
2. **Check firewall rules** that might block external connections
3. **Verify the VM has internet access**
4. **Check your organization's network policies**

## Storage and Data

### What storage options are supported?

**Current limitation**: Only Inductiva storage is currently supported in the alpha version.

**Coming soon**: GCP storage integration will be available in future releases.

### How do I handle large datasets?

Since only Inductiva storage is currently supported, you'll need to:
1. Upload your data to Inductiva storage first
2. Use the data in your simulations
3. Download results back to Inductiva storage

## CLI Usage

### How do I launch a BYOC machine from the command line?

```bash
inductiva task-runner launch byoc-gcp-machine \
    --provider gcp \
    --machine-type c2d-highcpu-8 \
    --spot
```

### How do I configure disk size via CLI?

```bash
inductiva task-runner launch byoc-gcp-machine \
    --provider gcp \
    --machine-type c2d-highcpu-8 \
    --data-disk-gb 100 \
    --spot
```

### How do I set custom auto-termination time via CLI?

```bash
inductiva task-runner launch short-lived-machine \
    --provider gcp \
    --machine-type c2d-highcpu-8 \
    --max-idle-time 5 \
    --spot
```

## Limitations and Future Features

### What are the current limitations?

**Alpha version limitations**:
- Single machine per group (no `num_machines` support)
- Cost reporting always shows 0
- Only Inductiva storage supported (GCP storage integration coming soon)

### What features are coming in future releases?

- **GCP storage integration** for direct access to your GCP storage buckets
- **Multi-machine scaling** within single machine groups
- **Improved cost reporting** and billing integration
- **Additional cloud providers** beyond GCP

### Is BYOC production-ready?

BYOC is currently in **alpha release**. While it's functional and secure, it has the limitations mentioned above. We recommend using it for development and testing, and monitoring the release notes for production-ready features.

```{banner_small}
:origin: launch_machines_on_gcp_faq
```
