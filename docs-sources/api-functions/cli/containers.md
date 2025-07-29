# inductiva **containers** [\[subcommands\]](#subcommands) [\[flags\]](#flags)

The `inductiva containers` command provides utilities for managing custom user containers. It allows users to convert Docker images into Apptainer-compatible `.sif` files and upload them to their Inductiva private storage for use in simulations.

##  Subcommands
### `convert`

Transforms a local or remote Docker image into a `.sif` file using Apptainer. 
This can be useful for users that want to see the conversion result and test the SIF file before uploading it to the remote storage. 

```bash
inductiva containers convert <IMAGE> <OUTPUT>
```
- `<IMAGE>` is the Docker image reference. Can be a Docher Hub reference (docker:// URL), a local image (name or ID) or a `.tar` archive exported from Docker.
- `<OUTPUT>` is the local path to save the resulting `.sif` file.

Sample usage:

```bash
# Convert a local image to SIF
inductiva containers convert my-image:latest ./my-image.sif

# Convert a Docker Hub image to SIF
inductiva containers convert docker://python:3.11-slim ./python.sif
```

### `upload`

Converts a Docker image to a `.sif` file and uploads it to the user's Inductiva private storage.

```bash
inductiva containers upload <IMAGE> [<OUTPUT_PATH>]
```
- `<IMAGE>` is the Docker image reference. Can be a Docher Hub reference (`docker://` URL), a local image (name or ID) or a `.tar` archive exported from Docker.
- `<OUTPUT_PATH>` is the path where the `.sif` file will be stored in Inductiva storage.
Defaults to: `my-containers/<image-name>.sif` if omitted.

Sample usage:

```bash
# Convert and upload a local Docker image
inductiva containers upload my-simulation-image

# Convert and upload a Docker Hub CFD image (SU2) with a custom storage path
inductiva containers upload docker://su2code/su2:latest my-containers/su2-cfd.sif
```

### `list (ls)` [\[flags\]](#flags-for-list)
List all container files stored in remote storage. If **no folder is provided**, it defaults to the `my-containers` directory.

```bash
inductiva containers list [<FOLDER>]
```

List container files from a specific folder:
```sh
inductiva containers list my-containers/project-a
```

Sample output:

```bash
$ inductiva containers ls
 NAME             SIZE        CREATION TIME
 container1.sif   200.00 MB   26/03, 16:41:14
 container2.sif   100.00 MB   26/03, 16:41:14
 container3.sif   300.00 MB   26/03, 16:41:14

Total storage size used:
        Volume: XXX GB
        Cost: YYY US$/month
```

<h4 id="flags-for-list">Flags</h4>

**`--max-results, -m`** (default:10)

Limits the number of results returned.

### `remove (rm)` [\[flags\]](#flags-for-remove)

Delete a container file from your Inductiva remote storage.
This action is **permanent** and cannot be undone — use with caution.

If no folder is provided, it defaults to `my-containers`.

```bash
inductiva containers remove [<FOLDER>]
```

Sample usage:

```bash
# Remove a container with confirmation
inductiva containers rm -n nginx.sif

# Remove a container without confirmation prompt
inductiva containers rm -n nginx.sif -y

# Remove a container from a specific folder
inductiva containers rm my-custom-folder -n my-container.sif
```

<h4 id="flags-for-remove">Flags</h4>

**`--name, -n`** (required)

The name of the container file to be removed.

---

**`--yes, -y`**

Skip the confirmation prompt and delete the container immediately.

## Flags
### `-h, --help`

Show help message and exit.

## Need Help?
Run the following command for more details:

```sh
inductiva containers --help
```

```{banner_small}
:origin: cli-containers
```
