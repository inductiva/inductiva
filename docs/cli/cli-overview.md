# CLI Overview and Functionalities

The Inductiva Command Line Interface (CLI) offers a streamlined and efficient way 
to interact with our platform directly from your terminal. Designed for ease of use, 
the CLI enables you to manage simulations, resources and logs in real-time, even 
while your simulations are running remotely.

## CLI Installation

The power of the Inductiva CLI lies in its simplicity, it's just a command away. 
Simply run on your terminal:

```bash
$ inductiva
```

## CLI Usage Guide

The Inductiva CLI streamlines the management of computational resources, enabling 
you to efficiently prepare and launch the necessary environments for your simulations. 

After installing you can choose from a list of available subcommands to run via 
the Inductiva (CLI) to manage your projects and interact with our API:

````{eval-rst}
.. tabs::

   .. tab:: `resources`

      `List and manage your computational resources. <./managing-resources.md>`_

   .. tab:: `tasks`

      Create and track your simulation tasks.

   .. tab:: `logs`

      Access real-time logs for ongoing tasks.

   .. tab:: `storage`

      Get an overview of your remote storage.

````

For detailed information on any subcommand, you can always run the `--help` or `-h` flag 
on your terminal to further understand its usage:

```bash
$ inductiva --help
```
