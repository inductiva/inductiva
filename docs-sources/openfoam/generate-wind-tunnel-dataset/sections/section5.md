## Postprocessing with Inductiva

Once all your OpenFOAM simulations have completed, the next step is often to extract useful physical quantities for analysis. Common postprocessing tasks include:

- Extracting **force coefficients** (drag, lift, moment, etc.)
- Generating **pressure maps** and **velocity fields**
- Visualizing **streamlines**, **vorticity**, or **flow separation**

In this section, we demonstrate how to automate one such task: **extracting force coefficients**, specifically, moment, drag, lift, front lift, and rear lift—from each simulation and compiling them into a single CSV file. This allows you to compare results across configurations, such as different wind speeds.

By running postprocessing as a cloud task with Inductiva:
- You **avoid downloading large datasets** to your local machine
- You **reuse existing simulation outputs** efficiently via `remote_assets`
- You maintain a **clean and reproducible workflow**
- You can **scale** to hundreds or thousands of simulations seamlessly

---

### Reusing Outputs with `remote_assets`

To make the process efficient, we don’t need to manually collect or download simulation results. Instead, we launch a new Inductiva task that directly reuses the output of the simulation tasks via the `remote_assets` feature.

Each previous task's output folder is referenced remotely in the cloud and accessed by the postprocessing script without needing to duplicate or transfer files.

> 🔗 For more information, check out: [Reusing Files from Other Tasks](https://inductiva.ai/guides/scale-up/reuse-files/reuse-files)

---

### Running the Postprocessing Task

We launch a lightweight Python task on a cloud machine and run a script (`postprocess.py`) that iterates through the output folders of all simulation tasks.

Below is the launcher script that:
- Collects the output files from the original simulation tasks
- Runs the postprocessing script in a `python:slim` container
- Waits for completion and downloads the results


```python
import inductiva

# Allocate cloud machine on Google Cloud Platform
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-highcpu-2",
    spot=True)

project = inductiva.projects.Project("openfoam_dataset")

remote_assets = [task.info.storage_output_path for task in project.get_tasks()]

python_image = inductiva.simulators.CustomImage(container_image="python:slim")

task = python_image.run(
    input_dir="input_files/",
    commands=["python postprocess.py"],
    on=cloud_machine,
    remote_assets=remote_assets)

task.wait()
cloud_machine.terminate()

task.download_outputs()
task.print_summary()
```

---

### The Postprocessing Script

The script scans each folder, reads the last line of `forceCoeffs.dat`, and extracts the relevant aerodynamic coefficients. It then compiles the results into a structured CSV file.


```python
import os
import csv

# Labels for coefficient data
coefficients_labels = ["Moment", "Drag", "Lift", "Front Lift", "Rear Lift"]
output_csv = "all_coefficients.csv"
rows = []

cwd = os.getcwd()

for folder in os.listdir(cwd):
    folder_path = os.path.join(cwd, folder)
    if not os.path.isdir(folder_path):
        continue

    coeffs_path = os.path.join(
        folder_path, "postProcessing", "forceCoeffs1", "0", "forceCoeffs.dat"
    )

    try:
        with open(coeffs_path, "r") as f:
            lines = [line.strip() for line in f if line.strip()]
            last_line = lines[-1]
            parts = last_line.split()

        values = list(map(float, parts[1:6]))
        row = {"Task id": folder}
        row.update(dict(zip(coefficients_labels, values)))
        rows.append(row)

    except Exception as e:
        print(f"Error in '{folder}': {e}")

# Write to CSV
if rows:
    fieldnames = ["Task id"] + coefficients_labels
    with open(output_csv, mode="w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
```

---

### Output

When the task completes, it produces an `all_coefficients.csv` file that looks like this:

```
Task id, Moment, Drag, Lift, Front Lift, Rear Lift
3b9n21xqqt97nbec07yzr6wzr, -0.034, 0.28, -1.34, -0.52, -0.82
d5r521g6igus8wh9c3yy8fbry, -0.056, 0.36, -1.65, -0.73, -0.92
...
```

This consolidated CSV is ready for visualization, comparison, statistical analysis, or even training a surrogate model.

---

### 🔍 Why Use a Separate Task?

Using a separate Inductiva task for postprocessing offers several advantages:

- ✅ **Reproducibility** – The script is versioned and can be rerun reliably
- ✅ **Portability** – Easy to adapt to new datasets or workflows
- ✅ **Efficiency** – No need to download or re-upload heavy simulation results
- ✅ **Scalability** – Supports processing large datasets without extra effort

Inductiva can simplify research by making high-performance computing more accessible and cost-effective.

```{banner_small}
:origin: openfoam
```
