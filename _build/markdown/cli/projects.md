# projects

## inductiva projects - CLI interface

```default
inductiva projects [-h] {download,list} ...
```

Projects management utilities.

The `inductiva projects` command offers tools to manage your projects  and access related information, such as listing all projects and  downloading their associated task output files.

### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

---

### inductiva projects download

```default
inductiva projects download [-h] [--output-dir OUTPUT_DIR] [--files FILES [FILES ...] |
                            --std]
                            project_name
```

The `inductiva projects download` command allows you to download the task outputs associated with a specific project. You can specify the files to download with the `--files` flag  or download the standard output and error files with the `--std` flag. If no files are specified, all files are downloaded.

#### Positional Arguments

* [**`project_name`**]() - Name of the project to download. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`--output-dir`**]() `OUTPUT_DIR` - Directory to save the downloaded files. (default: `None`)
* [**`--files`**]() `FILES` - Downloads the specified files from every task in the project. (default: `None`)
* [**`--std`**]() - Downloads the standard output and error files.

#### Examples

```bash

# Download `stderr.txt` and `system_metrics.csv` from all tasks in
# `my-swash-project`, and save them to the current directory.
$ inductiva projects download my-swash-project \
        --files stderr.txt system_metrics.csv \
        --output-dir .
Downloading simulation files to 9yhsu7cqcqemgzesuzueso1hn/outputs...
Downloading simulation files to e5ysb582dwofj3bw6233bae99/outputs...
Downloading simulation files to umpokj4e8rbx4379hsjdwfyuu/outputs...
Downloading simulation files to 1k3rvmf3e0cycsiiu62umg089/outputs...
Downloading simulation files to 9qj3fhe5tjt9bx8uvada9e8uy/outputs...
Partial download completed to 9yhsu7cqcqemgzesuzueso1hn/outputs.
Partial download completed to 1k3rvmf3e0cycsiiu62umg089/outputs.
Partial download completed to umpokj4e8rbx4379hsjdwfyuu/outputs.
Partial download completed to e5ysb582dwofj3bw6233bae99/outputs.
Partial download completed to 9qj3fhe5tjt9bx8uvada9e8uy/outputs.
Downloads completed.
```

---

### inductiva projects list (ls)

```default
inductiva projects list [-h] [-w [WATCH]]
```

The `inductiva projects list` command provides an overview of your projects. It lists all your available projects.

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-w`**]() `WATCH`, [**`--watch`**]() `WATCH` - Prompt the command every N seconds. (default: `None`)

#### Examples

```bash

$ inductiva projects list

NAME                                   NR_TASKS
openfoam-dambreak-2d-standard            5
openfoam-dambreak-3d-highres            22
xbeach-storm-surge-portugal-coast       10
openfast-wind-turbine                    4
xbeach-sediment-transport-validation     8
```

---
::docsbannersmall
::
