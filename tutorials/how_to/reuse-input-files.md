# Reuse Input Files Across Multiple Simulations

Running multiple simulations often involves reusing large input files, 
such as the bathymetry of a coastal area or object geometries for computational 
fluid dynamics (CFD). 

Instead of uploading these large files repeatedly for every simulation, 
Inductiva now offers a way to upload input files once and reuse them 
across multiple tasks. This not only saves time but also reduces costs.

In this tutorial, we’ll guide you through:
1. [Uploading input files to a remote directory.](#step-1-upload-input-files-to-remote-storage)
2. [Using these files in simulations.](#step-2-use-the-uploaded-files-in-simulations)
3. [Managing and maintaining remote files.](#step-3-maintain-and-manage-remote-files)

We've also included an [FAQ section](#faqs) to address common questions 
and help you get started quickly.

## Step 1: Upload Input Files to Remote Storage

You can upload files from either a local directory or a remote URL.

### From a Local Directory

Use the `upload` method to upload files or directories from your local 
system. Set the `local_path` parameter to the file or directory path you 
want to upload:

```python
inductiva.storage.upload(local_path="gromacs-input-example/",
    remote_dir="my_remote_directory",
)
```

### From a Remote URL

Use the `upload_from_url` method to upload files directly from a remote 
location. Set the `url` parameter to the remote file’s location:

```python
inductiva.storage.upload_from_url(url="https://storage.googleapis.com/inductiva-api-demo-files/test_assets/files.zip",
    remote_dir="my_remote_directory",
)
```
The `remote_dir` parameter specifies where the files will be stored remotely.

## Step 2: Use the Uploaded Files in Simulations

Once your files are uploaded, you can reference them in your simulations using 
the `remote_assets` parameter in the `simulator.run` method.

**Example with Gromacs Simulator**

```python
inductiva.storage.upload_from_url(url="https://storage.googleapis.com/inductiva-api-demo-files/test_assets/files.zip",
    remote_dir="my_remote_directory",
)
```

**Key Details:**

- The `remote_assets` parameter specifies the remote storage location 
where the input files are stored. This must match one of the directories 
you set as `remote_dir` in [Step 1](#step-1-upload-input-files-to-remote-storage).
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

machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

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
    remote_assets=["gromacs_bucket/file1.txt", "gromacs_bucket/file2.txt"])
```

## Step 3: Maintain and Manage Remote Files

You can list, clean, or manage your remote files directly through the 
Inductiva API or CLI.

### List Remote Files
- In CLI: `inductiva storage ls`
- With Python: `inductiva.storage.listdir()`

### Remove Files or Directories

- Remove an entire directory: `inductiva.storage.remove_workspace(remote_dir="gromacs_bucket")`
- Remove a single file from a remote directory: `inductiva.storage.remove_workspace(remote_dir="gromacs_bucket/file1.txt")`

### FAQs

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

*No, you don’t have to re-upload everything. Simply remove the specific file using `inductiva.storage.remove_workspace()` and upload the updated version. Alternatively, you can overwrite the file by uploading it via `input_dir`. Refer to question 10.

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