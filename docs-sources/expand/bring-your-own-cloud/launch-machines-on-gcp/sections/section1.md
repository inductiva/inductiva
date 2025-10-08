# How It Works

When you create a BYOC machine group, Inductiva creates a VM instance in your GCP account that runs a **task-runner container**. This container acts as a bridge between your GCP infrastructure and Inductiva's backend services. 

Here's the flow:

1. **VM Creation**: Inductiva creates a VM in your GCP project with the specified machine type and configuration
2. **Task-Runner Deployment**: The VM automatically launches a task-runner container that connects to Inductiva's backend
3. **Seamless Integration**: The task-runner handles communication with Inductiva's API and manages simulation execution
4. **Resource Management**: You maintain control over the VM lifecycle while Inductiva manages the simulation orchestration


> **Security safeguard**: Your GCP credentials and API keys are **never sent to Inductiva** and Inductiva **never has access to them**. All GCP operations (VM creation, deletion, metadata management) are performed entirely on your local machine using your authenticated gcloud CLI, ensuring your credentials remain secure and private.

> **What Inductiva Can See**: Even though we don’t have access to your GCP credentials, the task-runner running on your VM acts as a manager and reports machine information to Inductiva for pricing and monitoring purposes. This includes machine specifications (vCPUs, RAM, disk size), usage metrics (CPU usage, RAM usage, disk usage), runtime information (start time, last seen), live logs, and live file tracking for running simulations. This data helps Inductiva provide cost estimates and monitor machine health, but does not require access to your GCP credentials.


For more detailed information about the task-runner, see the [Task-Runner Guide](../../use-local-task-runner/index.md).

```{banner_small}
:origin: launch_machines_on_gcp_sec1b
```
