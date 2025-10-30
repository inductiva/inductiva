# containers

## inductiva containers - CLI interface

```default
inductiva containers [-h] {convert,list,ls,remove,rm,upload} ...
```

Manage custom simulation containers.

The `inductiva containers` command provides utilities for managing user-defined containers, including converting Docker images to Apptainer-compatible  `.sif` files and uploading them to your private Inductiva remote storage for use in simulations.

### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

---

### inductiva containers convert

```default
inductiva containers convert [-h] image output
```

The `inductiva containers convert` command converts a Docker image into  a Singularity Image Format (SIF, `.sif`) file using the Apptainer.

The Docker image can be specified as a Docker Hub URL, a local Docker image reference, or a `.tar` archive.

This can be useful for users that want to see the conversion result and  test the SIF file before uploading it to the remote storage.

#### Positional Arguments

* [**`image`**]() - Docker image to convert. Accepts a:
  	- local image name or ID (e.g., my-image:latest)
  	- Docker Hub reference URL (e.g., docker://username/image:tag)
  	- .tar archive exported from Docker (default: `None`)
* [**`output`**]() - The local path to save the resulting .sif file. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit

#### Examples

```bash

# Convert a local image to SIF
$ inductiva containers convert my-image:latest ./my-image.sif

# Convert a Docker Hub image to SIF
$ inductiva containers convert docker://python:3.11-slim ./python.sif
```

---

### inductiva containers list (ls)

```default
inductiva containers list [-h] [--max-results MAX_RESULTS] [folder]
```

The `inductiva containers list` command lists all container files in remote storage under the default containers folder (`my-containers/`), including their size and estimated cost.

You can also specify other container folders to list their contents.

#### Positional Arguments

* [**`folder`**]() - Path to a containers folder in remote storage. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`--max-results`**]() `MAX_RESULTS`, [**`-m`**]() `MAX_RESULTS` - Maximum number of container files to list. (default: `10`)

#### Examples

```bash

$ inductiva containers list
NAME             SIZE        CREATION TIME
container1.sif   200.00 MB   26/03, 16:41:14
container2.sif   100.00 MB   26/03, 16:41:14
container3.sif   300.00 MB   26/03, 16:41:14

Total storage size used:
        Volume: 0.59 GB
        Cost: 0.02 US$/month
```

---

### inductiva containers remove (rm)

```default
inductiva containers remove [-h] -n NAME [-y] [folder]
```

The `inductiva containers remove` command removes a specific container file from your Inductiva remote storage. If no folder is specified, defaults to `my-containers/`.

Use the flag `--yes` to skip confirmation prompts.

This action is irreversible and should be used with caution.

#### Positional Arguments

* [**`folder`**]() - Path to folder in remote storage. (default: `my-containers/`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-n`**]() `NAME`, [**`--name`**]() `NAME` - Name of the container file to remove. (default: `None`)
* [**`-y`**](), [**`--yes`**]() - Skip confirmation prompt and delete immediately.

#### Examples

```bash

# Remove a container with confirmation
$ inductiva containers rm -n nginx.sif

# Remove a container without confirmation prompt
$ inductiva containers rm -n nginx.sif -y

# Remove a container from a specific folder
$ inductiva containers rm my-custom-folder -n my-container.sif
```

---

### inductiva containers upload

```default
inductiva containers upload [-h] [-f] image [output_path]
```

Converts a Docker image (from Docker Hub, a local image, or a .tar file) into a SIF file using Apptainer, stores it in a temporary folder, and uploads that folder to your Inductiva remote storage, making it available for use with the Inductiva API.

#### Positional Arguments

* [**`image`**]() - Docker image reference. Accepts a:
  	- local image name or ID (e.g., python:3.11-slim)
  	- Docker Hub reference URL (e.g., docker://nginx:latest)
  	- .tar archive exported from Docker (default: `None`)
* [**`output_path`**]() - Optional output path for the .sif file in Inductiva remote storage (e.g., my-containers/nginx.sif). If omitted, defaults to: my-containers/<image-name>.sif. (default: `None`)

#### Options

* [**`-h`**](), [**`--help`**]() - show this help message and exit
* [**`-f`**](), [**`--overwrite`**]() - Overwrites the file in remote storage without asking.

#### Examples

```bash

# Convert and upload a local Docker image
$ inductiva containers upload my-simulation-image

# Convert and upload a Docker Hub CFD image (SU2) with a custom storage path
$ inductiva containers upload docker://su2code/su2:latest my-containers/su2cfd.sif
```

---
::docsbannersmall
::
