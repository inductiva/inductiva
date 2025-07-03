# containers

The inductiva containers command provides utilities for managing custom user containers. It allows users to convert Docker images into Apptainer-compatible .sif files and upload them to their Inductiva private storage for use in simulations.

# Usage

```bash
inductiva containers [-h] {convert,upload}
```

### Description
This command enables users to:

- Convert a local or remote Docker image into a Singularity Image Format (SIF) file.

- Converts a docker image to SIF file and uploads to Inductiva’s remote storage, making it available for use with the Inductiva API.

- List the containers default folder for visualization of container sizes and costs

###  Available Subcommands

- `convert` → Convert a Docker image into a `.sif` file.

- `upload` → Convert and upload a Docker image directly to your Inductiva storage.

- `list` → Displays the `my-containers` default folder.


## convert

The convert subcommand transforms a Docker image into a `.sif` file using Apptainer. 
This can be useful for users that want to see the conversion result and test the SIF file before uploading it to the remote storage. 

```bash
inductiva containers convert [-h] <image> <output>
```

### Positional Arguments
`image` → The Docker image to convert. Accepts:

- A local image name or ID (e.g., my-image:latest)

- A Docker Hub reference (e.g., docker://username/image:tag)

- A .tar archive exported from Docker.

`output` → The local path where the .sif file should be saved.


### Example

```bash
# Convert a local image to SIF
inductiva containers convert my-image:latest ./my-image.sif

# Convert a Docker Hub image to SIF
inductiva containers convert docker://python:3.11-slim ./python.sif
```

## upload

The upload subcommand both converts a Docker image to a .sif file and uploads it to the user's Inductiva private storage.


### Usage

```bash
inductiva containers upload [-h] <image> [output_path]
```

### Positional Arguments
`image` → The Docker image to convert. Accepts:

- A local image name or ID

- A Docker Hub reference (e.g., docker://nginx:latest)

- A .tar archive exported from Docker.


`output_path` (optional) → Path where the .sif file will be stored in Inductiva storage.
Defaults to: my-containers/<image-name>.sif if omitted.

### Example

```bash
# Convert and upload a local Docker image
inductiva containers upload my-simulation-image

# Convert and upload a Docker Hub CFD image (SU2) with a custom storage path
inductiva containers upload docker://su2code/su2:latest my-containers/su2-cfd.sif
```
## list

The list, or ls, subcommand will display the user's default container folder (my-containers).
Where all of the existing containers can be visualized and their size / cost to keep are displayed.

### Usage

```bash
inductiva containers list [-h]
```

### Positional Arguments

`folder` (optional) → Path to show the content.
`m` (optional) → max results to display.

### Example

```bash
➜  inductiva containers ls
 NAME             SIZE        CREATION TIME
 container1.sif   200.00 MB   26/03, 16:41:14
 container2.sif   100.00 MB   26/03, 16:41:14
 container3.sif   300.00 MB   26/03, 16:41:14

Total storage size used:
        Volume: XXX GB
        Cost: YYY US$/month
```

## remove

The remove (or rm) subcommand deletes a container file from your Inductiva remote storage. This action is irreversible and should be used with caution.


```bash
inductiva containers remove [-h] -n <name> [folder] [-y]
```


### Positional Arguments

`folder` (optional) → Folder path in remote storage. Defaults to my-containers.
`-n, --name` (required) → The name of the container file to remove.
`-y, --yes` (optional) → Skip confirmation prompt before deletion.

### Example

```bash
# Remove a container with confirmation
inductiva containers rm -n nginx.sif

# Remove a container without confirmation prompt
inductiva containers rm -n nginx.sif -y

# Remove a container from a specific folder
inductiva containers rm my-custom-folder -n my-container.sif
```

If the container is found, the CLI will prompt you for confirmation unless the --yes flag is used.

