# Inductiva: a Python package for scaling simulations on the Cloud

Welcome to the official Python library for the Inductiva API version 0.8
The Inductiva API allows running a set of open-source physical
simulators on the cloud, easily parallelizing simulations, each running
on hundreds of CPU cores.

Inductiva simplifies the complexities of cloud resource management, and software
configuration, offering a straightforward Python interface for running simulations
on state-of-the-art hardware. This allows scientists and engineers to focus their
time and energy on what matters: running simulations that solve real problems.

This documentation includes:
- A [Quick Recipes](./how_to/index.md) where we will guide you through
several common tasks and provides practical examples for you to try;
- A [list of the various open-source simulation packages](./simulators/overview.md)
that you can invoke via the API (if you have suggestions for adding more
simulation packages, please let us know by
[opening a GitHub issue](https://github.com/inductiva/inductiva/issues));
- A guide on [Inductiva’s Command Line Interface (CLI)](./cli/cli-overview.md), which
allows you to perform many tasks from your terminal, including listing available
computational resources and checking the status of tasks;
- A [User Reference](./api_reference/computational_resources/index.md) section 
that covers a wide variety of topics of interest, including information about
some key classes available in the API Client, a troubleshooting guide, information
about quotas and an [FAQ](./api_reference/faq.md).

If you have any questions or suggestions about the API please
[open an issue on the inductiva’s API Client GitHub repo](https://github.com/inductiva/inductiva/issues),
or contact us via [support@inductiva.ai](mailto:support@inductiva.ai).
