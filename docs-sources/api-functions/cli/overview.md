# Inductiva CLI

This section documents every command and flag available in Inductiva's command-line interface (CLI).

The Inductiva CLI helps you interact with the Inductiva platform directly from the terminal, enabling you to manage computational resources and tasks.

With the CLI, you can:
- Manage authentication and API keys
- Start, stop, list, and monitor [Tasks](../../how-it-works/tasks/index.md).
- Allocate and monitor [Computational Resources](../../how-it-works/machines/index.md).
- Manage projects and [storage files](../../how-it-works/storage/index.md)

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
    
    You are already logged in. Run `inductiva auth logout` if you want to log out. 
    Setting a new API Key will erase the existing one.
    To log in, you need an API Key. You can obtain it from your account at https://console.inductiva.ai/account.
Please paste your API Key here: 
```

## Need Help?

Use the `--help` or `-h` flag with **any command** for more details:

```sh
inductiva tasks --help
```
