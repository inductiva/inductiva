# 2D Ising Model with R and Inductiva

This tutorial shows how to run a parallel 2D ising-model Metropolis Monte Carlo simulation on R using the Inductiva API. 

<img src=_static/animation.gif></img>


In this tutorial we will: 

1. Prepare the R script
2. Setup the python script to call the Inductiva Python Client
3. Partition the temperatures across multiple machines and submit tasks
4. Collect and concatenate results
5. Plot magnetization curves


## Model description

The 2D Ising model is a simple yet powerful way to explore how local interactions give rise to large-scale order in magnetic materials. Imagine a square grid (lattice) of tiny “spins,” each of which can point either up (+1) or down (–1). Neighbors prefer to align—spins that point the same direction lower the system’s energy, while oppositely oriented neighbors increase it. The total energy of the lattice is simply the sum of all neighbor‐to‐neighbor interactions.

We drive the system forward using the Metropolis Monte Carlo algorithm: at each step we pick a random spin, calculate how flipping it would change the energy, and then decide probabilistically whether to accept that flip based on the temperature T. At low T, flips that increase energy are very unlikely, so the lattice quickly settles into a nearly uniform state (most spins aligned), producing a large overall magnetization. At high T, thermal agitation overwhelms the alignment preference and the lattice remains disordered, with roughly equal up and down spins.

Between these extremes lies a critical temperature $T_c$, where the system undergoes a phase transition: below T_c, it spontaneously “chooses” one of the two aligned states (positive or negative magnetization), while above T_c it stays essentially unmagnetized. Near T_c, the simulation exhibits rich behavior—very large clusters of aligned spins form and dissolve, and the time to reach equilibrium grows dramatically. In this tutorial, by sweeping T across a range of values and tracking the average magnetization over Monte Carlo sweeps, we will clearly see how order emerges and disappears in this classic model of statistical physics.

## R script

The `run_ising_model.R` is a simple command-line tool that simulates how a grid of magnetic “spins” flips and settles over time at different temperatures. You can pass in a few options when you run it:

1. **Pick your grid size**  
   Use `--N 100` to make a 100×100 grid of spins. Bigger grids give prettier patterns but take longer to compute.

2. **Choose your temperatures**  
   Use `--temps "1.5,2.0,2.5"` to run three temperatures in one go.  
   - Lower values (e.g. 1.5) mean “colder” spins that like to line up.  
   - Higher values (e.g. 2.5) add more random flipping.

3. **Decide how long to run**  
   - `--nsweeps 20000` makes the script attempt to flip every spin 20 000 times.  
   - `--burnin 2000` skips the first 2 000 sweeps so you only record data after the system has “warmed up.”  
   - `--thin 20` means “after burn-in, keep one measurement every 20 sweeps” to keep the output file manageable.

4. **What happens under the hood**  
   - Start with a random mix of up/down spins.  
   - Repeatedly pick each spin and ask:  
     1. “If I flip this spin, will its neighbors be happier (lower energy)?”  
     2. “If yes, always flip. If not, flip with a small chance that depends on **temperature**.”  
   - Over many sweeps, this ruleset naturally drives the grid toward full alignment at low T, or disorder at high T.

5. **What it records**  
   Every time you hit a “keep” sweep (after burn-in, every `thin` sweeps), it measures:  
   - **Magnetization**: how aligned the whole grid is (a number between –1 and +1).  
   - **Energy**: a single number that decreases as more spins agree with their neighbors.  

These measurements are written to a CSV with one row per snapshot:

`Temperature, SweepNumber, Magnetization, Energy`

> You can get the R script from this [link](bring-your-own-software/_static/run_ising_model.R).

## Python script

Below is a concise, step-by-step breakdown of how we drive the R simulation on Inductiva and retrieve the results.

---

### 1. Setup and configuration

We begin by importing the Inductiva client, defining where our R script lives (`input_dir`), and where to stash results (`results_dir`). We also list all the temperatures we want to try and set our simulation parameters.

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


## Visualization

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
