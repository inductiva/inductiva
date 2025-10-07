# Launching Machines on Your Own GCP Account

Inductiva's Bring Your Own Cloud (BYOC) feature allows you to launch and manage compute machines directly on your own Google Cloud Platform (GCP) account. This gives you full control over your compute resources while leveraging Inductiva's task management and simulation capabilities. Your GCP credentials never leave your local machine, ensuring maximum security and privacy.

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

**Security Note**: Your GCP credentials and API keys are **never sent to Inductiva** and Inductiva **never has access to them**. All GCP operations (VM creation, deletion, metadata management) are performed entirely on your local machine using your authenticated gcloud CLI, ensuring your credentials remain secure and private.

For more detailed information about the task-runner, see the [Task-Runner Guide](../use-local-task-runner/index.md).

## Prerequisites

Before you begin, ensure you have:

- **A Google Cloud Platform (GCP) account** with a valid billing account attached
- **Inductiva Python client** installed: See the [Installation Guide](../../../how-it-works/get-started/install-guide.md) for detailed instructions.

## Get Started

**[Installation and Setup](sections/section1b.md):** Installation steps, permissions, GCP CLI configuration, testing, and disclaimers.

**[Python Client Usage](sections/section2.md):** Create and manage BYOC machine groups using the Inductiva Python client.

**[Command Line Usage](sections/section3.md):** Launch and manage BYOC machine groups using the Inductiva CLI.

**[Troubleshooting and Limitations](sections/section4.md):** Common issues, solutions, and current limitations of the alpha release.

## Disclaimer

When using BYOC, you are responsible for:
- **All costs** incurred by VMs running in your GCP account
- **Managing and monitoring** your compute resources
- **Ensuring compliance** with your organization's policies
- **Properly terminating** machines to avoid unexpected charges

**Note:** Inductiva provides auto-termination (default: 3 minutes of idle time) as a safety feature to help prevent forgotten VMs. Inductiva provides support on a best-effort basis and cannot guarantee resolution of all issues.

**Inductiva is not responsible for:**
- Unexpected charges from VMs left running in your account
- Data loss or security incidents on your GCP infrastructure
- Compliance issues with your organization's policies
- Any damage to your GCP resources or data

**Recommendation:** Always monitor your GCP console and set up billing alerts to track usage and costs.

```{toctree}
:hidden:
sections/section1b.md
sections/section2.md
sections/section3.md
sections/section4.md
```

```{banner}
:origin: launch_machines_on_gcp
```
