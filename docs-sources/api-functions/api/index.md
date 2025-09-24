# Python Client

The **Inductiva Python Client** is a **library** that transforms the Inductiva API requests into simple Python code. Instead of dealing with complex API endpoints manually, we provide Python functions and classes that let you manage computational resources, run simulations, organize projects, and retrieve results — all from a simple Python script.

With the Python Client you can focus on your simulation design and analysis rather than API mechanics, making it easy to integrate high-performance computational simulations into your existing Python scripts!

## Classes & Methods Overview

The table below provides an overview of the Python Client's main classes and their available methods. Most of these methods have corresponding CLI commands that offer the same functionality through the command line interface.

For a more comprehensive understanding on when to use each interface and how they work together, see our [Interfaces with the API](http://inductiva.ai/guides/how-it-works/building-blocks/index) guide.

| Class        | Methods                                 | CLI                                             | Resource Guide                                                   |
|----------------------|---------------------------------------------|-----------------------------------------------------------|------------------------------------------------------------------|
| [`benchmarks`](inductiva.benchmarks)               | `add_run`, `add_task`, `download_outpus`, `export`, `get_tasks`, `run`, `set_default`, `terminate`, `wait`                           | --                              | [Benchmark Guide](https://inductiva.ai/guides/scale-up/benchmark/index)        |
| [`projects`](inductiva.projects)               | `add_task`, `download_outputs`, `get_tasks`, `wait`                           | [Projects CLI](../cli/projects.md)                              | [Projects Guide](https://inductiva.ai/guides/scale-up/projects/index)        |
| [`resources`](inductiva.resources)              | `active_machines_to_str`, `can_start_resource`, `estimate_cloud_cost`, `get_available_mpi_slots`, `get_mpi_config`, `start`, `terminate`, `get`, `get_by_name`, `estimate_machine_cost`, `get_available_machine_types`        | [Resources CLI](../cli/resources.md)                              | [Resources Guide](https://inductiva.ai/guides/how-it-works/machines/index)                |
| [`storage`](inductiva.storage)            | `copy`, `download`, `export`, `export_to_aws_s3`, `get_signed_urls`, `get_space_used`, `get_zip_contents`, `listdir`, `multipart_upload`, `remove`, `upload`, `upoload_from_url`                | [Storage CLI](../cli/storage.md)                          | [Storage Guide](https://inductiva.ai/guides/how-it-works/intro/data_flow)            |
| [`tasks`](inductiva.tasks)          | `close_stream`, `download_inputs`, `download_outputs`, `get_info`, `get_input_url`, `get_machine_type`, `get_metadata`, `get_output_info`, `get_output_url`, `get_position_in_queue`, `get_simulator_name`, `get_status`, `get_storage_path`, `get_total_time`, `is_failed`, `is_running`, `is_terminal`, `kill`, `last_modified_files`, `list_files`, `print_summary`, `remove_remote_files`, `set_metadata`, `tail_files`, `wait`, `wait_for_status`      | [Tasks CLI](../cli/tasks.md)                      | [Tasks Guide](https://inductiva.ai/guides/how-it-works/tasks/index)          |
| [`templating`](inductiva.templating)               | `render_dir`              | --      | [Templating Guide](https://inductiva.ai/guides/scale-up/parallel-simulations/templating)  |
| [`users`](inductiva.users)         | `get_costs`, `get_info`, `get_quotas`                  | [Users CLI](../cli/user.md)  | -- |

---

## Set Up & Authentication

To get started, follow our [Installation Guide](https://inductiva.ai/guides/how-it-works/get-started/install-guide) to install the [Inductiva Python package](https://pypi.org/project/inductiva/) and set up your environment.

```{banner}
:origin: api
```