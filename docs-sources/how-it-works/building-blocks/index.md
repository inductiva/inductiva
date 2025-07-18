# Building Blocks

Inductiva provides an API to easily run large-scale computational simulations on the cloud.

To cater to different workflows and user preferences, we offer three ways to interact with our API: the **Python Client**, the **Command-Line Interface (CLI)**, and the **Web Console**. While each tool can be used independently, they are designed to work together seamlessly, providing a flexible and unified experience for managing your simulations.

This page explains the purpose of each component and how they fit together to form the Inductiva ecosystem.

## How It All Connects

At the heart of Inductiva is the API that manages all your cloud resources, simulation tasks, and storage. The Python Client, CLI, and Web Console are simply different interfaces that communicate with these resources. This means any action you take in one tool will be immediately reflected in the others.


````{eval-rst}
.. mermaid::
   :align: center

   graph TD
      subgraph User Interfaces
         A[Python Client]
         B[CLI]
         C[Web Console]
      end

      subgraph "Inductiva Platform"
         API(Inductiva API)
      end

      subgraph "Cloud Infrastructure"
         D[Computing Resources]
         E[Cloud Storage]
      end

      A -- Programmatic Access --> API
      B -- Terminal Commands --> API
      C -- Visual Interaction --> API

      API -- Manages --> D
      API -- Manages --> E

      style A fill:#f9f,stroke:#333,stroke-width:2px
      style B fill:#ccf,stroke:#333,stroke-width:2px
      style C fill:#9cf,stroke:#333,stroke-width:2px

      linkStyle 0,1,2 stroke-width:2px,fill:none,stroke:orange;
      linkStyle 3,4 stroke-width:2px,fill:none,stroke:green;

````


The following diagram illustrates this relationship:

This page explains the purpose of each component and how they fit together to form the Inductiva ecosystem.

Inductiva provides three complementary ways to interact with our API, each designed for different use cases and user preferences. Whether you're a engineer building complex simulation workflows or a scientist stydying complex phenomenon

Welcome to the official Python library for the Inductiva API version 0.10
The Inductiva API allows running a set of open-source physical
simulators on the cloud, easily parallelizing simulations, each running
on hundreds of CPU cores.

Inductiva simplifies the complexities of cloud resource management, and software
configuration, offering a straightforward Python interface for running simulations
on state-of-the-art hardware. This allows scientists and engineers to focus their
time and energy on what matters: running simulations that solve real problems.

This documentation includes:

- An [introduction to our API](../intro/how_it_works.md), highlighting its functionalities and key features.
- A guide on [Inductiva’s Command Line Interface (CLI)](../CLI/index.md), which
allows you to perform many tasks from your terminal, including listing available
computational resources and checking the status of tasks;
- A [User Reference](./computational_resources/index.md) section
that covers a wide variety of topics of interest, including information about
some key classes available in the API Client, a troubleshooting guide, information
about quotas and an [FAQ](../FAQ/index.md).

If you have any questions or suggestions about the API please
[open an issue on the inductiva’s API Client GitHub repo](https://github.com/inductiva/inductiva/issues),
or contact us via [support@inductiva.ai](mailto:support@inductiva.ai).
