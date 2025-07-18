# Inductiva Python Client

The Inductiva **Python Client** provides a programmatic interface to the Inductiva API, enabling you to manage computational resources, run simulations, organize projects, and retrieve results — all from within a **Python script**.

Our Python Client is a library that wraps the Inductiva API calls into Python functions and classes. With this, you can launch simulations, manage resources, and process results using straightforward Python code — no need to manage complex API endpoints manually!

## Python Client Classes Overview

| Class        | Methods                                 | CLI                                             | Resource Guide                                                   |
|----------------------|---------------------------------------------|-----------------------------------------------------------|------------------------------------------------------------------|
| [`benchmarks`](inductiva.benchmarks)               | `add_run`, `add_task`, `download_outpus`, `export`, `get_tasks`, `run`, `set_default`, `terminate`, `wait`                           | --                              | [Benchmark Guide](../../scale-up/benchmark/index.md)        |
| [`projects`]()               | `add_task`, `download_outputs`, `get_tasks`, `wait`                           | [Projects CLI](../cli/projects.md)                              | --        |
| [`resources`]()              | `active_machines_to_str`, `can_start_resource`, `estimate_cloud_cost`, `get_available_mpi_slots`, `get_mpi_config`, `start`, `terminate`, `get`, `get_by_name`, `estimate_machine_cost`, `get_available_machine_types`        | [Resources CLI](../cli/resources.md)                              | [Resources Guide](../../how-it-works/machines/index.md)                |
| [`storage`]()            | `copy`, `download`, `export`, `export_to_aws_s3`, `get_signed_urls`, `get_space_used`, `get_zip_contents`, `listdir`, `multipart_upload`, `remove`, `upload`, `upoload_from_url`                | [Storage CLI](../cli/storage.md)                          | [Storage Guide](../../how-it-works/cloud-storage/index.md)            |
| [`tasks`]()          | `close_stream`, `download_inputs`, `download_outputs`, `get_info`, `get_input_url`, `get_machine_type`, `get_metadata`, `get_output_info`, `get_output_url`, `get_position_in_queue`, `get_simulator_name`, `get_status`, `get_storage_path`, `get_total_time`, `is_failed`, `is_running`, `is_terminal`, `kill`, `last_modified_files`, `list_files`, `print_summary`, `remove_remote_files`, `set_metadata`, `tail_files`, `wait`, `wait_for_status`      | [Tasks CLI](../cli/tasks.md)                      | [Tasks Guide](../../how-it-works/tasks/index.md)          |
| [`templating`]()               | `render_dir`              | --      | [Templating Guide](../../scale-up/parallel-simulations/templating.md)  |
| [`users`]()         | `get_costs`, `get_info`, `get_quotas`                  | [Users CLI](../cli/user.md)  | -- |

---

## Set Up & Authentication

To get started, follow our [Installation Guide](https://inductiva.ai/guides/how-it-works/get-started/install-guide) to install the [Inductiva Python package](https://pypi.org/project/inductiva/) and set up your environment.

```{toctree}
---
caption: Python Client
maxdepth: 2
hidden: true
---
benchmarks <inductiva.benchmarks>
projects <inductiva.projects>
resources <inductiva.resources>
storage <inductiva.storage>
simulators <inductiva.simulators>
tasks <inductiva.tasks>
templating <inductiva.templating>
users <inductiva.users>
```
