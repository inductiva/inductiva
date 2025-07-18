# Inductiva CLI

The **Inductiva CLI** helps you interact with the Inductiva platform directly from the terminal, enabling you to manage computational resources and tasks.

This section documents every command and flag available in Inductiva's command-line interface (CLI).

With the CLI, you can:
- Manage authentication and API keys
- Start, stop, list, and monitor [Tasks](../../how-it-works/tasks/index.md).
- Allocate and monitor [Computational Resources](../../how-it-works/machines/index.md).
- Manage projects and [storage files](../../how-it-works/storage/index.md)

## CLI Command Overview

| Command Group        | Subcommands                                 | API Reference                                             | Resource Guide                                                   |
|----------------------|---------------------------------------------|-----------------------------------------------------------|------------------------------------------------------------------|
| `auth`               | `login`, `logout`                           | --                              | [Authentication Guide](../../how-it-works/get-started/install-guide.md)        |
| `tasks`              | `download`, `info`, `kill`, `last-modified-file`, `list`, `list-files`, `tail`, `top`        | [Tasks API](https://inductiva.ai/guides/api-functions/api/inductiva.tasks)                              | [Tasks Guide](../../how-it-works/tasks/index.md)                |
| `task-runner`          | `launch`, `remove`      | --                     | [BYOH Guide](../../expand/use-local-task-runner/index.md)          |
| `projects`           | `list`, `download`        | [Projects API](https://inductiva.ai/guides/api-functions/api/inductiva.projects)                        | --       |
| `storage`            | `list`, `download`, `export`, `remove`, `size`                | [Storage API](https://inductiva.ai/guides/api-functions/api/inductiva.storage)                          | [Storage Guide](../../how-it-works/cloud-storage/index.md)            |
| `resources`          | `start`, `cost`, `info`, `available`, `list`, `terminate`      | [Resources API](https://inductiva.ai/guides/api-functions/api/inductiva.resources)                      | [Machines Guide](../../how-it-works/machines/index.md)          |
| `simulators`               | `list`              | [Simulators API](https://inductiva.ai/guides/api-functions/api/inductiva.simulators)      | [Simulators Guide](../../how-it-works/simulators/index.md)  |
| `containers`         | `list`, `upload`, `remove`                  | --  | [BYOH Guide](../../expand/bring-your-own-software/index.md) |

---

### Notes

- Use `--help` with any command for more detailed usage and available flags.
- `containers` operates on top of the remote storage system, designed specifically for `.sif` container files.
- For safety, commands like `resources terminate` and `containers remove` include a `--yes` flag to bypass confirmation prompts.


---

## Set Up & Authentication

The CLI is automatically installed when you install
Inductiva's [Python Client](../api/index.md). However, before using the CLI,
you need to authenticate with your Inductiva account.

```sh
inductiva auth login
```

When prompted, copy paste your personal API key that is availalbe
from your [User Account](<https://console.inductiva.ai/account/profile>)
page in the Web Console.

```sh
     ___  _   _  ____   _   _   ____  _____  ___ __     __ _
    |_ _|| \ | ||  _ \ | | | | / ___||_   _||_ _|\ \   / // \
     | | |  \| || | | || | | || |      | |   | |  \ \ / // _ \
     | | | |\  || |_| || |_| || |___   | |   | |   \ V // ___ \
    |___||_| \_||____/  \___/  \____|  |_|  |___|   \_//_/   \_\
    
    To log in, you need an API Key. You can obtain it from your account at https://console.inductiva.ai/account.
    Please paste your API Key here: 
```

## Need Help?

Use the `--help` or `-h` flag with **any command** for more details:

```sh
inductiva tasks --help
```
