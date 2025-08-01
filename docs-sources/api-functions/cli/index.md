# Inductiva CLI

The **Inductiva CLI** helps you interact with the Inductiva API directly from the terminal, enabling you to manage computational resources, tasks, and much more.

This section documents every command and flag available in Inductiva's command-line interface (CLI).

## CLI Command Overview

The table below shows the available CLI commands alongside their corresponding Python Client APIs and documentation guide. Most functionality is available through both interfaces, allowing you to choose the tool that best fits your workflow. The CLI is ideal for quick operations, while the Python Client offers a more programmatic control and integration with your code.

For detailed guidance on when to use each interface and how they work together, see our [Interfaces with the API](http://inductiva.ai/guides/how-it-works/building-blocks/index) guide.

| Command        | Python Client                                             | Resource Guide                                                   |
|----------------------|-----------------------------------------------------------|------------------------------------------------------------------|
| [`auth`](auth.rst)               | --                              | [Authentication Guide](https://inductiva.ai/guides/how-it-works/get-started/install-guide)        |
| [`user`](user.rst)               | [Users API](https://inductiva.ai/guides/api-functions/api/inductiva.users)                              | --        |
| [`tasks`](tasks.rst)              | [Tasks API](https://inductiva.ai/guides/api-functions/api/inductiva.tasks)                              | [Tasks Guide](https://inductiva.ai/guides/how-it-works/tasks/index)                |
| [`task-runner`](task-runner.rst)          | --                     | [BYOH Guide](https://inductiva.ai/guides/expand/use-local-task-runner/index)          |
| [`projects`](projects.rst)           | [Projects API](https://inductiva.ai/guides/api-functions/api/inductiva.projects)                        | [Projects Guide](https://inductiva.ai/guides/scale-up/projects/index)       |
| [`storage`](storage.rst)            | [Storage API](https://inductiva.ai/guides/api-functions/api/inductiva.storage)                          | [Storage Guide](https://inductiva.ai/guides/how-it-works/intro/data_flow)            |
| [`resources`](resources.rst)          | [Resources API](https://inductiva.ai/guides/api-functions/api/inductiva.resources)                      | [Resouces Guide](https://inductiva.ai/guides/how-it-works/machines/index)          |
| [`simulators`](simulators.rst)               | [Simulators API](https://inductiva.ai/guides/api-functions/api/inductiva.simulators)      | [Simulators Guide](../../how-it-works/simulators/index.rst)  |
| [`containers`](containers.rst)         | --  | [BYOS Guide](https://inductiva.ai/guides/expand/bring-your-own-software/index) |

### Notes

- Use `--help` with any command for more detailed usage and available flags.
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

```{banner}
:origin: cli
```