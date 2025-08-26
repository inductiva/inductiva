# Visualizing Results with Smokeview

[Smokeview](https://github.com/firemodels/smv) is a tool for visualizing the results of Fire Dynamics Simulator (FDS) simulations.
In this tutorial, we’ll show how to run Smokeview **directly on Inductiva** after your FDS simulation finishes, using the "Town House Kitchen Fire" example from the [Smokeview repository](https://github.com/firemodels/smv/tree/SMV-6.9.1/Verification/Visualization).

## Prerequisites

Before starting, you should have already [run your first FDS simulation](setup-test.md) using the Inductiva API.
You’ll also need the input files required for this tutorial: `thouse5.fds` and `thouse5_movies.ssf`, with minor modifications.

You can download the ready-to-use files [here](https://storage.googleapis.com/inductiva-api-demo-files/fds-tutorials/SmokeviewExampleKitchenFire.zip).

## Running the Simulation with Visualization

To use Smokeview in **scripted mode** via the Inductiva API:

* Include a `.ssf` file (Smokeview script) in your input directory.
* Optionally include a `.ini` file to define visualization preferences.
* Set the `smokeview_case_name` argument to match the base name of your simulation (i.e., the `.smv` file produced by FDS).

Smokeview will be automatically run using:

```bash
smokeview -runscript <smokeview_case_name>
```

### Generating Videos from Smokeview Rendered Frames

Smokeview can render sequences of PNG images using the `RENDERALL` command in your `.ssf` file.
These image sequences can then be automatically compiled into videos using the `smokeview_video_sources` argument.

For example:

If your `.ssf` file contains:

```
RENDERDIR .
...
RENDERALL
 1
thouse5_tslice
...
RENDERALL
 1
thouse5_smoke3d
```

Smokeview will generate files like:

```
thouse5_tslice_0001.png, thouse5_tslice_0002.png, ...
thouse5_smoke3d_0001.png, thouse5_smoke3d_0002.png, ...
```

To generate videos from these, pass the prefixes to `smokeview_video_sources`:

```python
smokeview_video_sources=["thouse5_tslice", "thouse5_smoke3d"]
```

---

## Full Example: Running the "Town House Kitchen Fire"

```python
import inductiva

# Allocate a cloud machine
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2d-standard-2",
    spot=True)

# Initialize the FDS simulator
fds = inductiva.simulators.FDS(version="6.9.1")

# Run the simulation with visualization enabled
task = fds.run(
    input_dir="/Path/to/SmokeviewKitchenFireExample",
    sim_config_filename="thouse5.fds",
    smokeview_case_name="thouse5",
    smokeview_video_sources=["thouse5_tslice", "thouse5_smoke3d"],
    on=cloud_machine)

# Wait for simulation to complete and download results
task.wait()
cloud_machine.terminate()
task.download_outputs()
task.print_summary()
```

> **Note**: Set `input_dir` to the local path containing your `.fds`, `.ssf`, and any related files.

## Output

The following two MP4 videos will be produced:

<p align="center">
  <img src="_static/thouse5_smoke3d.gif" alt="Smoke3D visualization" width="45%" style="display:inline-block; margin-right: 10px;">
  <img src="_static/thouse5_tslice.gif" alt="T-slice visualization" width="45%" style="display:inline-block;">
</p>

## Key Takeways

You can use Inductiva to streamline your FDS post-processing with Smokeview, by using the following arguments:
- `smokeview_case_name`: Must match the base name of the `.smv` and `.ssf` files (without extension).
- `smokeview_video_sources`: PNG filename prefixes used to generate videos for each rendered sequence.

