# Building Blocks

Inductiva provides an API to easily run large-scale simulations on the cloud.

To accomodate different workflows and user preferences, we offer three ways to interact with our API: the [Python Client](https://inductiva.ai/guides/api-functions/api/index), the [Command-Line Interface (CLI)](https://inductiva.ai/guides/api-functions/cli/overview), and the [Web Console](https://console.inductiva.ai/dashboard). While each tool can be used independently, they are designed to work together, providing a unified experience for managing your simulations.

This page explains the purpose of each component and how they all fit together.

## Architecture Overview

The Inductiva API centrally manages all your cloud resources, simulation tasks, and storage. The Python Client, CLI, and Web Console are the different interfaces you can use to communicate with these resources.

The following diagram illustrates this relationship:

![Building Blocks](../_static/building_blocks.png)

## The Three Building Blocks

### Python Client: _a library that transforms API requests into simple code_

The **Python Client** is a library to control the Inductiva API programmatically within a Python script.

**Best For:**
- Automating complex or repetitive workflows 
- Programmatically manage resources and data storage.

```{seealso}
Check the complete Python Client documentation `here <https://inductiva.ai/guides/api-functions/api/index>`_.
```

### Command Line Interface: _for quick terminal operations_

The **Inductiva CLI** provides a fast and efficient way to interact with the API directly from the terminal, without needing to write a full script.

**Best For:**
- Check the status of a task or list all running machines.
- Download outputs or view logs for a specific task.
- Manually start or terminate a resource.

```{seealso}
Check the complete CLI documentation `here <https://inductiva.ai/guides/api-functions/cli/overview>`_.
```

### Web Console: _an intuitive graphical interface for visualization and management_

The **Web Console** is a graphical interface to visually monitor tasks and resources, analyze costs, see statistics and analytics, and manage your account. It requires **no programming or command-line knowledge**, making it the most suited tool for a high-level overview.

**Best For:**
- Visually monitor resources.
- Analyze costs and manage your credits and account settings.
- Perform urgent actions like killing a task or terminating a machine group with a few clicks.

## How They Fit Together

Inductiva's ecosystem is designed for maximum flexibility. With all three interfaces built on the same API, you can seamlessly move between Python scripts, terminal commands, and the web interface as your workflow evolves.

Consider a typical workflow:

1.  **Launch with Python:** define and launch simulations using a **Python script**.
2.  **Monitor with the CLI:** while the script is running, use the **CLI** in your terminal to quickly check the status of your tasks (`inductiva tasks list`) or view the live logs (`inductiva tasks tail <TASK_ID>`).
3.  **Visualize in the Console:** At any given moment, you can open the [Web Console](https://console.inductiva.ai/dashboard) to visually inspect the simulation outputs and metrics (cost, time, system metrics, for example).



| | AUTOMATE & SCRIPT | LAUNCH SIMULATIONS | MONITOR & OPERATE | VISUAL ANALYTICS | MANAGE ACCOUNT |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Python Client** | ✅ | ✅ | ✅ | ❌ | ❌ |
| **CLI** | ❌ | ✅ | ✅ | ❌ | ❌ |
| **Web Console** | ❌ | ❌ | ✅ | ✅ | ✅ |