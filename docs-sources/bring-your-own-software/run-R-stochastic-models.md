# How to Run Statistical Ensembles with Inductiva

In this recipe-style tutorial, you’ll learn how to run *statistical ensembles* in parallel using Inductiva, so you can turn any parameter sweep or distribution into dozens (or hundreds!) of simultaneous simulations. We’ll demonstrate with a classic 2D Ising-model Metropolis Monte Carlo in R, but the same pattern works for *any* function or simulation script that takes parameters and spits out results.

<img src=_static/animation.gif></img>

> **Why Inductiva?**  
> Instead of looping through parameters one at a time, Inductiva lets you launch a group of machines, run tasks concurrently, and collect outputs automatically. This means you can get ensemble results in **minutes** instead of **hours**.


**What you’ll do in this tutorial:**

1. **Package your simulation**  
   Convert your simulation code (in this case, an R script) to accept parameters from the command line and output results to CSV.

2. **Define your ensemble**  
   Specify a list or distribution of parameters (e.g., temperatures, seeds, reaction rates) to explore.

3. **Launch an Inductiva machine group**  
   Allocate as many workers as ensemble tasks.

4. **Dispatch tasks in parallel**  
   Send each parameter set as an individual job. Inductiva manages queuing, scaling, retries, and logging.

5. **Visualize ensemble results**  
   Plot summary statistics like means, variances, or distributions using a few lines of code.



> **Model link:** If you need the 2D Ising-model details, see [2D Ising model description »](https://en.wikipedia.org/wiki/Ising_model)  


## Prepare your ensemble script

First, wrap your ensemble script logic into a command-line script that:

- **Reads** one or more parameter values from the CLI  
- **Runs** your script over those parameters  
- **Writes** its results to an output (in this case a csv)

Below is an example layout for [run_ising_model.R](bring-your-own-software/_static/run_ising_model.R). In this tutorial it takes a list of “temperature” values, but you can swap in *any* numeric parameter (reaction rate, random seed, boundary condition,).

```bash
Rscript run_ising_model.R \
  --temps "1.5,2.0,2.5" \
  --grid-size 100 \
  --n-sweeps 20000 \
  --burn-in 2000 \
  --thin 20 \
  --output results.csv
```

| Flag           | Description                                                      |
| -------------- | ---------------------------------------------------------------- |
| `--temps`      | Comma-separated values to sweep (e.g. `"1.5,2.0,2.5"`).          |
| `--grid-size`  | Size of your simulation grid (e.g. `100` ⇒ 100×100).             |
| `--n-sweeps`   | Total sweeps per run (e.g. `20000`).                             |
| `--burn-in`    | Number of initial sweeps to skip before recording (e.g. `2000`). |
| `--thin`       | Record one snapshot every *n* sweeps after burn-in (e.g. `20`).  |
| `--output`     | CSV filename for your results (e.g. `results.csv`).              |


>💡 Pro tip: Any script that reads --param-list and writes a CSV with columns ParamValue, Sweep, Output1, Output2, … can slot into this same Inductiva workflow.

## Launching your ensemble with Inductiva

Suppose we want to examine how minor perturbations around a few base temperatures affect the behavior of our Ising model. With Inductiva, we can sample a normal distribution around each base temperature, run each simulation in parallel, and collect the full set of results—without manually managing infrastructure.

Here's how to do that in a single, self-contained Python script:


---

### 1. Setup and configuration

Import the necessary libraries, point to your R script folder (input_dir), and pick where to save results (results_dir). Define the base temperatures, the number of samples around each, and the Normal‐sampling spread (sigma).

```python
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import inductiva

INPUT_DIR = "./input"
RESULTS_DIR = "./results"
OUTPUT_ROOT = "inductiva_output"
PROJECT_NAME = "temperature_ensemble"

BASE_TEMPERATURES = [1.5, 2.0, 2.5]
NUM_SAMPLES = 5
TEMPERATURE_SIGMA = 0.2
RNG_SEED = 123

NSWEEPS = 10000
BURNIN = 1000
N_SPINS = 100
THIN = 10

# Ensure results directory exists
os.makedirs(RESULTS_DIR, exist_ok=True)
```

### 2. Sample a Normal cloud around each T0

For each base temperature T₀, draw n_samples points from 𝒩(T₀, σ²). That gives you one “central” run at T₀ plus a small spread around it.

```python
np.random.seed(RNG_SEED)
jobs_by_base = {
    T0: sorted(np.random.normal(loc=T0, 
                                scale=TEMPERATURE_SIGMA, 
                                size=NUM_SAMPLES))
         for T0 in BASE_TEMPERATURES
}
```

### 3. Spin up the cloud machines

Compute the total number of runs (one for each T₀ plus its samples) and use Inductiva API to spin up that many machines. Here we use an ElasticMachineGroup, which grows automatically based on how many tasks are being executed, you can control how many machines you want to have using the `max_machines` parameter.

```python
total_jobs = sum(1 + len(samples) for samples in jobs_by_base.values())

machine_group = inductiva.resources.ElasticMachineGroup(
    machine_type="c2-standard-8",
    spot=True,
    min_machines=1,
    max_machines=total_jobs,
)

# Use the R docker container simulator
sim = inductiva.simulators.CustomImage(container_image="docker://r-base:latest")

# Define a project for the current tasks
project = inductiva.projects.Project(project_name)
```

### 4. Submit every temperature as its own task

Loop over each base T₀ and its sampled neighbors, round each temperature, build the Rscript command, and submit. We can attach metadata so we can later pull out exactly which runs belong to which T₀.


```python
for T0, samples in jobs_by_base.items():
    # include the base temperature itself plus the sampled ones
    for temp in [T0] + samples:
        temp_rounded = round(temp, 3)
        output_filename = f"ising_base{T0:.1f}_sample{temp_rounded:.3f}.csv"

        command = (f'Rscript run_ising_model.R '
                   f'--N {N_SPINS} '
                   f'--temps "{temp_rounded}" '
                   f'--nsweeps {NSWEEPS} '
                   f'--burnin {BURNIN} '
                   f'--thin {THIN} '
                   f'--output {output_filename}')

        task = simulator.run(
            input_dir=INPUT_DIR,
            on=machine_group,
            commands=[command],
            project=PROJECT_NAME,
        )

        task.set_metadata({
            "base_temperature": str(T0),
            "sample_temperature": str(temp_rounded),
            "output_file": output_filename,
            "task_id": task.id,
        })

        print(
            f"Submitted task: T0={T0} temp={temp_rounded} → {output_filename}")
```

### 5. Wait for completion and download

Block until every task in the project finishes, then pull down all their CSV outputs into the default inductiva_output/... folder.
In this case, the output folder wil be `inductiva_output/<PROJECT_NAME>/`.


```python
# Wait for all tasks in the project to end
project.wait()

# Download all outputs from the project
project.download_outputs()
```

At this point we'll have all of the CSV per temperature, ready to merge and visualize.


## Visualize ensemble results

Finally, for each base temperature T₀, collect its CSVs, merge them, and plot the mean magnetization ⟨M⟩ vs sweep:

```python
# Gather all metadata from the project
all_tasks = project.get_tasks()
metadata_list = [task.get_metadata() for task in all_tasks]

for T0 in BASE_TEMPERATURES:

    # 1) Find all CSVs for this T0
    csv_paths = []
    for meta in metadata_list:
        if meta["base_temperature"] != str(T0):
            continue
        tid = meta["task_id"]
        fname = meta["output_file"]
        path = os.path.join(OUTPUT_ROOT, PROJECT_NAME, tid, "outputs", fname)
        csv_paths.append(path)

    # 2) Read & concatenate
    df = pd.concat((pd.read_csv(p) for p in csv_paths), ignore_index=True)

    # 3) Plot ⟨M⟩ vs sweep for each temperature
    plt.figure()
    for temp in sorted(df["T"].unique()):
        mean_M = df[df["T"] == temp].groupby("sweep")["M"].mean()
        plt.plot(mean_M.index, mean_M.values, label=f"T={temp}")

    plt.xlabel("Sweep")
    plt.ylabel("Magnetization ⟨M⟩")
    plt.title(f"Ising Magnetization for T₀={T0}")
    plt.legend()
    plt.tight_layout()

    plot_file = os.path.join(RESULTS_DIR, f"magnetization_T0_{T0}.png")
    plt.savefig(plot_file)
    plt.close()

    print(f"✔ Saved plot for T₀={T0} → {plot_file}")
```

<img src=_static/combined_magnetization.png></img>


And that completes our end-to-end workflow:

1. Prepared the R Monte Carlo script
2. Orchestrated parallel runs via Python + Inductiva
3. Collected every run's output
4. Merged and visualized the results

We can now observe how magnetization evolves across sweeps for each base temperature and its nearby samples.
