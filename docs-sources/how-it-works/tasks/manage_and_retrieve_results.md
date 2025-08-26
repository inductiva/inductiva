# Manage Task Results

The **output** of a simulation is its most important artifact. The Inductiva API provides tools for managing where these results are stored remotely and for downloading them to your local machine. Remember that every simulation you run creates a `Task` object, which is the key to accessing and managing these results.

This guide covers two main operations:
 * Controlling the remote directory where simulation outputs are saved.
 * Downloading results to your local machine, either in full or partially.

For a general introduction to `Task` objects, please refer to the [main Tasks documentation](index.md).

## Controlling Remote Storage Location

By default, all outputs from a simulation are saved in your remote storage space in a directory named after the unique `task.id`.

For example, running a simulation will generate a `Task` and store its results in a folder named with that task's ID:

```python
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2-standard-4")
cloud_machine.start()

simulator = inductiva.simulators.REEF3D()

task = simulator.run(
    input_dir="path-to-directory-with-input-files-for-reef3d",
    on=cloud_machine)
print(task.id)  # will print i4ir3kvv62odsfrhko4y8w2an

# Terminate the machine group
cloud_machine.terminate()
```

You can verify this using the Inductiva CLI to list the contents of your remote storage. The output will show a directory named `i4ir3kvv62odsfrhko4y8w2an/`:

```bash
$ inductiva storage list --max-results 2 --order-by creation_time --sort-order desc

       NAME                             SIZE           CREATION TIME
       i4ir3kvv62odsfrhko4y8w2an/       21.67 MB        07 Feb, 13:47:03

```

### Customizing the Remote Directory Name

To specify a custom remote directory name, use the `storage_dir` argument in the simulator's run method.

```python
task = simulator.run(
    input_dir="path-to-directory-with-input-files-for-reef3d"
    storage_dir="my_reef3d_simulation")
```

Listing your storage contents again will now show the new custom directory:

```bash
$ inductiva storage list --max-results 2 --order-by creation_time --sort-order desc

       NAME                             SIZE           CREATION TIME
       my_reef3d_simulation/            21.67 MB       07 Feb, 14:57:11
       i4ir3kvv62odsfrhko4y8w2an/       21.67 MB       07 Feb, 13:47:03

```

> Check out the complete remote storage documentation [here](https://inductiva.ai/guides/how-it-works/cloud-storage/cloud-storage)


## Downloading Results Locally
The primary way to download results in a Python script is with the `task.download_outputs()` method.

```python
task = inductiva.tasks.Task("i4ir3kvv62odsfrhko4y8w2an")
task.download_outputs()
```

````sh
Downloading simulation files to inductiva_output/i4ir3kvv62odsfrhko4y8w2an/outputs/output.zip...
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 56.5M/56.5M [00:01<00:00, 49.5MB/s]
Uncompressing the files to inductiva_output/i4ir3kvv62odsfrhko4y8w2an/outputs...
````

For a quick way to get results from the terminal, use the `inductiva tasks download` CLI command with the task ID.

By default, all output files will be download into a local folder structure `inductiva_output/<task_id>/`.

#### Customizing the Download Path

You can override the default download location in two ways:

1. For a specific download, pass the `output_dir` argument to `download_outputs()`: 

```python
task.download_outputs(output_dir="my-sim-results")
```

or using the CLI, pass the `--dir` argument:

```sh
inductiva tasks download i4ir3kvv62odsfrhko4y8w2an --dir my-sim-results
````

This will download the files to `my-sim-results/`.

2. For the entire session, use `inductiva.set_output_dir()` to change the global parent directory. Setting it to `None` will download directly to the current working directory.

```python
>>> import inductiva
# set the parent directory name
>>> inductiva.get_output_dir()
'inductiva_output'
>>> inductiva.set_output_dir("my_vault")
>>> inductiva.get_output_dir()
'my_vault'

# download the outputs again
>>> output_dir = task.download_outputs()
Downloading simulation outputs to my_vault/i4ir3kvv62odsfrhko4y8w2an/output.zip.
100%|██████████████████████████████████████| 5.04M/5.04M [00:00<00:00, 13.3MB/s]
Uncompressing the outputs to my_vault/i4ir3kvv62odsfrhko4y8w2an/
>>> print(ouptut_dir)
PosixPath('my_vault/i4ir3kvv62odsfrhko4y8w2an')

# Unset the parent directory name so that all downloads are done to the current
# working directory
>>> inductiva.set_output_dir(None)
>>> output_dir = task.download_outputs()
Downloading simulation outputs to i4ir3kvv62odsfrhko4y8w2an/output.zip.
100%|██████████████████████████████████████| 5.04M/5.04M [00:00<00:00, 13.3MB/s]
Uncompressing the outputs to i4ir3kvv62odsfrhko4y8w2an/
>>> print(ouptut_dir)
PosixPath('i4ir3kvv62odsfrhko4y8w2an')
```

### Downloading Specific Files

To download only a subset of the output files, use the `filenames` argument and provide a list of the files you need. This can save time and bandwidth.

````python
output_dir = task.download_outputs(filenames=["stdout.txt", "stderr.txt"]
                                       output_dir="my_outputs")
````

````sh
Downloading simulation outputs to my_outputs/output.zip.
100%|██████████████████████████████████████| 1.04k/1.04k [00:00<00:00, 671kB/s]
Uncompressing the outputs to my_outputs/
````

or using the CLI:

```sh
inductiva tasks download i4ir3kvv62odsfrhko4y8w2an --filenames stdout.txt stderr.txt --dir my_outputs
````

Checking the contents of the destination folder confirms that only the specified files were downloaded:

```bash
$ ls -1 my_outputs
stderr.txt
stdout.txt
```


