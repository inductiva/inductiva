# Reuse Files Across Multiple Simulations

Running multiple simulations often involves reusing large files, 
such as the bathymetry of a coastal area or object geometries for computational 
fluid dynamics (CFD). 

Instead of uploading these large files repeatedly for every simulation, 
Inductiva now offers a way to upload files once and reuse them 
across multiple tasks. This not only saves time but also reduces costs.

Additionally, you can use the outputs of one task as inputs for another, 
enabling easy checkpointing and reducing unnecessary data transfers.

In this tutorial, we’ll guide you through:
1. [Uploading input files to a remote directory.](#upload-input-files-to-remote-storage)
2. [Managing and maintaining remote files.](#maintain-and-manage-remote-files)
3. [Reusing task outputs in simulations.](#reuse-task-outputs-in-simulations)

We've also included an [FAQ section](#faqs) to address common questions 
and help you get started quickly.

## Upload Input Files to Remote Storage

You can upload files from either a local directory or a remote URL.

### From a Local Directory

Use the `upload` method to upload files or directories from your local 
system. Set the `local_path` parameter to the file or directory path you 
want to upload:

```python
inductiva.storage.upload(
    local_path="gromacs-input-example/",
    remote_dir="my_remote_directory")
```

### From a Remote URL

Use the `upload_from_url` method to upload files directly from a remote 
location. Set the `url` parameter to the remote file’s location:

```python
inductiva.storage.upload_from_url(
    url="https://storage.googleapis.com/inductiva-api-demo-files/test_assets/files.zip",
    remote_dir="my_remote_directory",
)
```
The `remote_dir` parameter specifies where the files will be stored remotely.

### Use the Uploaded Files in Simulations

Once your files are uploaded, you can reference them in your simulations using 
the `remote_assets` parameter in the `simulator.run` method.

**Example with Gromacs Simulator**

```python
inductiva.storage.upload_from_url(
    url="https://storage.googleapis.com/inductiva-api-demo-files/test_assets/files.zip",
    remote_dir="my_remote_directory",
)
```

**Key Details:**

- The `remote_assets` parameter specifies the remote storage location 
where the input files are stored. This must match one of the directories 
you set as `remote_dir` in [previous step](#upload-input-files-to-remote-storage).
- The `input_dir` parameter can still be used for local files. If no `remote_assets` 
are provided, the input files will be read from the local `input_dir`.
- If both `remote_assets` and `input_dir` are provided, and files with 
the same name exist in both locations, the files from `input_dir` will take priority.
- You only need to provide one of these parameters (`remote_assets` or `input_dir`).

### Use Multiple Remote Inputs

The `remote_assets` parameter accepts a list, allowing you to specify 
multiple remote files or directories:

```python
import inductiva

machine_group = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2-standard-4")

commands = [
    "gmx solvate -cs tip4p -box 2.3 -o conf.gro -p topol.top",
    ("gmx grompp -f energy_minimization.mdp -o min.tpr -pp min.top -po min.mdp "
     "-c conf.gro -p topol.top"),
    "gmx mdrun -s min.tpr -o min.trr -c min.gro -e min.edr -g min.log",
    ("gmx grompp -f positions_decorrelation.mdp -o decorr.tpr -pp decorr.top "
     "-po decorr.mdp -c min.gro"),
    ("gmx mdrun -s decorr.tpr -o decorr.trr -x  -c decorr.gro -e decorr.edr "
     "-g decorr.log"),
    ("gmx grompp -f simulation.mdp -o eql.tpr -pp eql.top -po eql.mdp "
     "-c decorr.gro"),
    ("gmx mdrun -s eql.tpr -o eql.trr -x trajectory.xtc -c eql.gro -e eql.edr "
     "-g eql.log"),
]

gromacs = inductiva.simulators.GROMACS()

task = gromacs.run(
    input_dir=None,
    commands=commands,
    on=machine_group,
    remote_assets=[
        "gromacs_bucket/file1.txt",
        "gromacs_bucket/file2.txt"])
```

```{banner_small}
:origin: scale_up_reuse_files_upload_remote_storage
```

## Maintain and Manage Remote Files

You can list, clean, or manage your remote files directly through the 
Inductiva API or CLI.

### List Remote Files
- In CLI: `inductiva storage ls`
- With Python: `inductiva.storage.listdir()`

### Remove Files or Directories

- Remove an entire directory: `inductiva.storage.remove(remote_path="gromacs_bucket")`
- Remove a single file from a remote directory: `inductiva.storage.remove(remote_path="gromacs_bucket/file1.txt")`

```{banner_small}
:origin: scale_up_reuse_files_manage_remote_files
```

## Reuse Task Outputs in Simulations

You can reuse the output of a previous task in a new simulation by including its `storage_output_path` in the `remote_assets` parameter. This is especially useful when chaining simulations, for example, using a minimized structure, topology file, or configuration generated by a prior task.


### What Happens Internally?

When you add a previous task’s `storage_output_path` to the `remote_assets` list:
- All files from that task’s output directory are made available in your new task.
- These files are automatically copied into your new task’s working directory under a subdirectory named after the original task’s ID.

For example, if the previous task has ID `2s9njoiy59jzr96m0sbbsmgwr`, then its files will be available in the folder: `2s9njoiy59jzr96m0sbbsmgwr/` within the new task’s workspace.

### Example 1: Reusing One Task’s Output

```python
previous_task = inductiva.tasks.Task("<task_id>")
task = gromacs.run(
    input_dir=None,
    commands=commands,
    on=machine_group,
    remote_assets=[previous_task.info.storage_output_path])
```

This makes the files from `<task_id>` available in a folder named `<task_id>/` inside the new task’s workspace.

### Example 2: Reusing Multiple Tasks

You can reuse the outputs of multiple previous tasks:

```python

task = gromacs.run(
    input_dir=None,
    commands=commands,
    on=machine_group,
    remote_assets=[
        previous_task_1.info.storage_output_path,
        previous_task_2.info.storage_output_path])
```

Each task’s output will be accessible in a folder named after that task’s ID:
```
<task_id_1>/...
<task_id_2>/...
```

### How to Reference Output Files

Let’s say the task stored in the `previous_task` variable produced a `topol.top` file, and you want to use it in your next simulation.

Because its output will be placed in a folder named after the task’s ID (i.e., `previous_task.id`) inside your new task’s workspace, your command should look like this:


```python
commands = [
    f"gmx solvate -cs tip4p -box 2.3 -o conf.gro -p {previous_task.id}/topol.top",
    ...
]
```
Make sure to use an f-string (`f"..."`) so that `previous_task.id` gets correctly inserted into the command string.


### How to Reference Task Outputs in a Config File
In some simulations, you may need to inject the path of a file from a previous task into a config file, rather than referencing it directly in a command. For example, in OpenFAST, to reuse a TurbSim-generated wind field, you must point to its .bts file in the `InflowWind.dat` file.

#### Step 1: Make the Config File a Template
Open the `InflowWind.dat` file and replace the hardcoded path with a template variable. For example:

Change this line:

```"./TurbSim"   FileName_BTS  - Name of the Full field wind file to use (.bts)```

to 

```{{bts_filepath}}   FileName_BTS  - Name of the Full field wind file to use (.bts)```

Then rename the file from: `InflowWind.dat` to `InflowWind.dat.jinja`.

This tells the platform to treat the file as a Jinja2 template that will be rendered with specific parameters.

#### Step 2: Inject the Path in Your Script
Assuming you previously ran TurbSim and saved its task object as `turbsim_task`, the wind field `TurbSim.bts` file will be located inside the task’s output.

You can now render the template like this:

```python
import inductiva

inductiva.TemplateManager.render_dir(..., bts_filepath=f"{turbsim_task.id}/TurbSim.bts")
```

> 💡 **Tip:** Want to learn more about using templates for simulation input files?
> Check out the [Inductiva Templating Guide](https://inductiva.ai/guides/scale-up/parallel-simulations/templating).

```{banner_small}
:origin: scale_up_reuse_outputs
```

## FAQs

**1. What file types can I upload?** 

*Any file type can be uploaded and reused.*

---

**2. Where are my files stored in the remote storage?**

*Your files are stored in your personal area within Inductiva’s filesystem.*

---

**3. Can I upload multiple files at once?**

*Yes, but only when you upload a local directory using `inductiva.storage.upload`.*

---

**4. What happens if I try to upload a file to a remote directory that already contains a file with the same name?**

*Inductiva will show an error. You need to remove the existing file before uploading the new one.*

---

**5. Can I use files from different remote directories in a single task?**

*Yes, you can specify multiple remote files and directories in the `remote_assets` parameter.*

---

**6. How do I confirm that my upload was successful?**

*You can use `inductiva.storage.listdir()` in Python or `inductiva storage ls` in the CLI to check the contents of your remote directory.*

---

**7. Can I upload a zip file?**

*Yes, you can upload zip files, but they won’t be automatically unzipped. We’re working on adding this feature in the future.*

---

**8. I need to change one file in a directory I uploaded. Do I need to re-upload the entire directory?**

*No, you don’t have to re-upload everything. Simply remove the specific file using `inductiva.storage.remove()` and upload the updated version. Alternatively, you can overwrite the file by uploading it via `input_dir`. Refer to question 10.

---

**9. Can I track the progress of my file upload?**

*Yes. If you’re uploading from your local system, a progress bar will appear. For remote uploads, you can use `inductiva.storage.listdir()` or `inductiva storage ls` to check progress.*

---

**10. What happens if files in `input_dir` and `remote_assets` have the same name?**

*Files in `input_dir` will take priority and overwrite files with the same name in `remote_assets`. This ensures that locally provided files always override remote files.*

*Example:*

```python
task = gromacs.run(
input_dir="local_folder/",
    commands=commands,
    on=machine,
    remote_assets=["remote_folder/"]
)
```
*If `object.obj` exists in both `local_folder` and `remote_folder`, the simulation will use the file from local_folder.*

---

**11. What happens if I remove or update a remote directory while it’s being used by a task?**

*Once the task starts, the remote assets are copied to the task. Any changes won’t affect the ongoing simulation.*

---

**12. Does using remote assets improve task performance?**

*Yes, reusing remote assets reduces the need to upload large files repeatedly, cutting down task startup time.*

---

**13. Can I use remote files directly in `remote_assets` without uploading them first?**

*Not yet, but we’re planning to add this feature soon.*


---

**14. Can I reuse a single file from another task output?**

*Yes, just add the `storage_path` of the task to the `remote_assets` list to reuse the output file.*

```{banner_small}
:origin: scale_up_faqs
```