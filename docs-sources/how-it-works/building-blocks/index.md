# Building Blocks

Inductiva provides an API to easily run large-scale computational simulations on the cloud.

To accomodate different workflows and user preferences, we offer three ways to interact with our API: the **[Python Client](../../api-functions/api/index.md)**, the **[Command-Line Interface (CLI)](../../api-functions/cli/overview.md)**, and the **Web Console**. While each tool can be used independently, they are designed to work together, providing a unified experience for managing your simulations.

This page explains the purpose of each component and how they all fit together.

## Architecture Overview

At the heart of Inductiva is the API that manages all your cloud resources, simulation tasks, and storage. The Python Client, CLI, and Web Console are the different interfaces we provide to communicate with these resources. This means any action you take in one tool will be immediately reflected in the others.

The following diagram illustrates this relationship:

![Building Blocks](../_static/building_blocks.png)

## The Three Building Blocks

### Python Client: _a library that transforms API requests into simple Python code_

The **Python Client** is a library that allows you to control the Inductiva API programmatically within a Python script.

**Best For:**
- Automate complex or repetitive workflows 
- Programmatically manage resources and data storage.
- Perform custom pre- and post-processing steps around your simulations.

### Command Line Interface: _for quick terminal operations_

The **Inductiva CLI** provides a fast and efficient way to interact with the API directly from the terminal. Suited for quick operations, allowing to manage resources and check on tasks without writing a full script.

**Best For:**
- Check the status of a task or list all running machines.
- Download outputs or view logs for a specific task.
- Manually start or terminate a resource.


### Web Console: _an intuitive graphical interface for visualization and management_

The **Web Console** is a graphical interface (UI) that provides a visual way to monitor tasks and resources, analyze costs, see statistics and analytics and manage your account. It requires **no programming or command-line knowledge**, making it the perfect tool for getting a high-level overview, performing quick actions, and managing billing.

**Best For:**
- Visually monitor running tasks and computational resources.
- Analyze costs and manage your credits and account settings.
- Perform urgent actions like killing a task or terminating a machine group with a few clicks.

## How They Fit Together: A Unified Experience

The true power of Inductiva comes from using these tools together. Because they all interact with the same API, you can mix and match them to suit your needs at any given moment.

Consider a typical workflow:

1.  **Launch with Python:** You define and launch simulations using a **Python script** because it allows for complex logic and templating.
2.  **Monitor with the CLI:** While the script is running, you can use the **CLI** in your terminal to quickly check the status of your tasks (`inductiva tasks list`) or view the live logs (`inductiva tasks tail <TASK_ID>`).
3.  **Visualize in the Console:** At any given moment, you can open the **Web Console** to visually inspect the graphical outputs, check the total cost of the run, and download the final results to your simulation tasks.


| TOOL | AUTOMATE & SCRIPT | LAUNCH SIMULATIONS | MONITOR & OPERATE | VISUAL ANALYTICS | MANAGE ACCOUNT |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Python Client** | ✅ | ✅ | ✅ | ❌ | ❌ |
| **CLI** | ❌ | ✅ | ✅ | ❌ | ❌ |
| **Web Console** | ❌ | ❌ | ✅ | ✅ | ✅ |