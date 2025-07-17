# inductiva **projects** [\[subcommands\]](#subcommands) [\[flags\]](#flags)
The `inductiva projects` command provides utilities to manage and
retrieve information about your projects.

Projects in Inductiva serve as folders for organizing tasks.
Using the CLI, you can list your projects and retrieve the
outputs of all tasks in a specific projects.

## Subcommands

### `list (ls)`
List all existing projects associated with your account.

```sh
inductiva projects list
```

or using the alias:

```sh
inductiva projects ls
```

### `download` [\[flags\]](#flags-for-download)
Download the task files of a specific project. By default, the outputs of  the project's tasks will be
download to a directory named `inductiva_output/{TASK_ID}/outputs`.

```sh
inductiva projects download <PROJECT_NAME>
```

<h4 id="flags-for-download">Flags</h4>

**`--output-dir=<directory>`**

Specify the directory to save the downloaded files.

---

**`--files=<filename>`**

Specify the files from every project's task you want to download.

For example, you can download the `stderr.txt` and `system_metrics.csv` of every task in the `default` project by doing:

```sh
$ inductiva projects download default --files stderr.txt system_metrics.csv
    Downloading simulation files to inductiva_output/g4z2navzczzcb66oqbkyfoqdn/outputs...
    Downloading simulation files to inductiva_output/0x6u6ccj3mmhzpa03xnxknrog/outputs...
    Partial download completed to inductiva_output/g4z2navzczzcb66oqbkyfoqdn/outputs.
    Partial download completed to inductiva_output/0x6u6ccj3mmhzpa03xnxknrog/outputs.
    ...
```

---

**`--std`**

Download the standard output and error files.

## Flags
### `-h, --help`

Show help message and exit.

## Need Help?
Run the following command for more details:

```sh
inductiva projects --help
```
