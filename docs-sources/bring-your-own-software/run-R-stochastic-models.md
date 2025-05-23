# How to Run Statistical Ensembles with Inductiva

In this recipe-style tutorial, you’ll learn how to run *statistical ensembles* in parallel using Inductiva, so you can turn any parameter sweep or distribution into dozens (or hundreds!) of simultaneous simulations. We’ll demonstrate with a classic 2D Ising-model Metropolis Monte Carlo in R, but the same pattern works for *any* function or simulation script that takes parameters and spits out results.

<img src=_static/animation.gif></img>

> **Why Inductiva?**  
> Instead of looping over parameters one by one, Inductiva lets you spin up a machine group, dispatch each task concurrently, and collect all outputs automatically, so you get your ensemble results back *minutes* instead of *hours*.

**What you’ll do in this tutorial:**

1. **Package your simulation**  
   Wrap your code (here, an R script) so it reads parameters from the command line and writes its results to CSV.

2. **Define your ensemble**  
   Specify a list or distribution of parameters (e.g. temperatures, initial seeds, reaction rates) that you want to explore.

3. **Launch an Inductiva machine group**  
   Spin up as many workers as tasks in your ensemble.

4. **Dispatch tasks in parallel**  
   Send each parameter set off as its own job; Inductiva handles queuing, scaling, retries, and logging.

5. **Gather and merge outputs**  
   Automatically download every CSV, stitch them together, and import into R (or Python) for analysis.

6. **Visualize ensemble results**  
   Plot summary statistics (means, variances, distributions) across your ensemble in a few lines of code.


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

Below is a concise, step-by-step breakdown of how you can drive the R simulation on Inductiva and retrieve the results.

---

### 1. Setup and configuration

We begin by importing the Inductiva client, defining where our R script lives (`input_dir`), and where to stash results (`results_dir`). We also list all the temperatures we want to try and set our ensemble parameters.

```python
import os, math, inductiva

# Local folder with run_ising_model.R
input_dir   = "./input"

# Temperatures to sweep and simulation settings
all_temps = [1.5, 1.8, 2.0, 2.2, 2.4, 2.8, 3.0, 3.2, 3.6, 4.0]
N, nsweeps, burnin, thin = 100, 10000, 1000, 10
```

### 2. Launch your machine group

Here we tell Inductiva Client to spin up as many machines as we have temperatures (one per T). These machines will each run their own R script concurrently.

```python
mg = inductiva.resources.ElasticMachineGroup(
    machine_type="c2d-highcpu-4",
    spot=True,
    min_machines=1,
    max_machines=len(all_temps),
)
```

### 3. Assign one temperature per machine

Since we have one machine for each temperature, we simply pair each T with a machine by slicing the list into single-element “chunks.”

```python
temp_chunks = [T for T in all_temps]
```

### 4. Submit the R jobs

For each temperature, we build the appropriate Rscript command and submit it. We track each returned task object alongside the expected output filename.


```python
sim = inductiva.simulators.CustomImage(container_image="docker://r-base:latest")
project = inductiva.projects.Project("2d-ising-model")

for idx, temps in enumerate(temp_chunks, start=1):
    out_csv = f"ising_results_chunk{idx}.csv"
    cmd = (
        f'Rscript run_ising_model.R '
        f'--N {N} --temps "{temps}" '
        f'--nsweeps {nsweeps} --burnin {burnin} --thin {thin} '
        f'--output {out_csv}'
    )
    task = sim.run(input_dir=input_dir, 
                   on=mg, 
                   commands=[cmd],
                   project=project.name,
    )
    print(f"→ Submitted T={temps} on task {task.id}")
```

### 5. Wait for completion and download

Finally, we block on each task until it finishes, then download its CSV into the default inductiva output folder.
In this case, the output folder wil be `inductiva_output/2d-ising-model/`.


```python
# Wait for all tasks in the project to end
project.wait()

# Download all outputs from the project
project.download_outputs()
```

At this point you’ll have all of the CSV per temperature, ready to merge and visualize. Each step is clear and keeps logic flow obvious: configure → launch machines → assign work → submit jobs → collect outputs.


## Visualize ensemble results

Now that all the tasks have finished, we can use their outputs to generate a visualization of the magnetization vs. sweep. Here’s how:

```python

# Read all .csv files from the output
pattern = os.path.join("inductiva_output", project.name, "*", "outputs",
                       "ising_results_chunk*.csv")
files = glob.glob(pattern)
csv_files = sorted(files)

# Read all of the csv output files
df = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)

plt.figure()
for T in sorted(df["T"].unique()):
    meanM = df[df["T"] == T].groupby("sweep").M.mean()
    plt.plot(meanM.index, meanM.values, label=f"T={T}")

plt.xlabel("Sweep")
plt.ylabel("Magnetization ⟨M⟩")
plt.title("2D Ising Magnetization vs. Sweep")
plt.legend()
plt.tight_layout()

plot_path = os.path.join(results_dir, "magnetization.png")
plt.savefig(plot_path)
```


<img src=_static/magnetization.png></img>


And with that, you have a complete end-to-end tutorial:

1. Prepared the R Monte Carlo script
2. Orchestrated parallel runs via Python + Inductiva
3. Collected every chunk’s output
4. Merged and visualized the results

You should now see how magnetization grows (or decays) over sweeps at each temperature.
