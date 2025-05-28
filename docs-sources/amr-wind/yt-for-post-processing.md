# AMR-Wind Post-Processing
Getting started with post-processing AMR-Wind simulation results is straightforward using [yt](https://yt-project.org/doc/visualizing/plots.html#slice-plots), a powerful Python package designed for analyzing and visualizing volumetric simulation data. It supports AMReX plotfiles out of the box and is ideal for scripting and automation.

Alternatively, if you prefer an interactive graphical interface, Paraview also offers native support for AMReX plotfiles and provides an easy way to explore and visualize your simulation results.

In this tutorial, we’ll guide you step-by-step through using yt to load your AMR-Wind simulation data, generate visualizations like slice plots, and automate the creation of animations to better understand your results.

## Using yt for Post-Processing
**yt** is ideal for scripting and automated analysis of AMR data. With just a few lines of Python, you can generate slice plots, projections, and derived fields quickly.

### Step 1: Install yt 
If you don’t have yt installed yet, you can add it easily with pip:

```bash
pip install yt
```

### Step 2: Load Data and Create a Slice Plot
Here’s a simple example to load a plotfile and generate a 2D slice of the x-velocity field along the z-axis:

```python
import yt

ds = yt.load("plt01000") # Load a time step of your choice from the provided path
slc = yt.SlicePlot(ds, "z", ("boxlib", "velocityx")) # Create a slice plot of x velocity normal to z-axis
slc.set_cmap(("boxlib", "velocityx"), "plasma") # Choose colormap
slc.set_log(("boxlib", "velocityx"), False) # Use linear scale (log scale is default)
slc.annotate_title("Velocity field slice")
slc.save("xvelocity_plot.png")
```

### 3. Automate Plotting Across Multiple Time Steps
You can easily loop over all `plt00000` files (time steps), extracting a data slice and saving it as a 2D plot.

```python
import glob, os
import yt

plotfiles = sorted(glob.glob("input/path_to_plt_files/plt*"))
for pf in plotfiles:
    ds = yt.load(pf)
    slc = yt.SlicePlot(ds, "z", ("boxlib", "mag_vorticity"))
    slc.set_cmap(("boxlib", "mag_vorticity"), "plasma")
    slc.set_log(("boxlib", "mag_vorticity"), False)
    slc.annotate_title(f"Time = {float(ds.current_time):.2f} s")
    slc.save() # Saves plot files in the same folder as the plt files

```

### Creating Animations from Plotfiles
To visualize temporal evolution, you can generate a series of slice plots and compile them into an animation (GIF). Below is a sample script demonstrating this process:

```python
import yt
import os
import glob
import imageio

def plot_variable_slices(data_path, output_path, variable, case_title):
    """
    Plots 2D z-slices of a given variable over time from AMR-Wind plotfiles using yt.

    Args:
        data_dir (str): Path to directory containing AMR-Wind plotfiles.
        output_subdir (str): Subdirectory (relative to data_dir) to save output images.
        variable (str): Variable to plot (e.g., 'velocityx', 'mag_vorticity').
        case_title (str): Case name or description to use in plot titles.
    """
    vmin, vmax = 0, 50  # Color scale limits

    os.makedirs(output_path, exist_ok=True)

    plotfiles = sorted(glob.glob(os.path.join(data_path, "plt*")))

    for pf in plotfiles:
        ds = yt.load(pf)

        # Create slice plot
        slice_plot = yt.SlicePlot(ds, "z", ("boxlib", variable))
        slice_plot.set_cmap(("boxlib", variable), "plasma")
        slice_plot.set_log(("boxlib", variable), False)
        slice_plot.set_zlim(("boxlib", variable), vmin, vmax)

        time_val = float(ds.current_time)
        slice_plot.annotate_title(f"{case_title}  |  Time = {time_val:.1f} s")
        slice_plot.hide_axes()

        # Save to file
        basename = os.path.basename(pf)
        slice_plot.save(os.path.join(output_path, f"{basename}.png"))

def create_gif_from_images(image_dir, output_gif, duration=0.5):
    """
    Creates a GIF from a series of PNG images.

    Args:
        image_dir (str): Directory containing PNG images.
        output_gif (str): Output path for the GIF.
        duration (float): Duration between frames in seconds.
    """
    image_paths = sorted(glob.glob(os.path.join(image_dir, "*.png")))
    images = [imageio.imread(img) for img in image_paths]
    imageio.mimsave(output_gif, images, duration=duration)

if __name__ == "__main__":
    # User parameters
    case_title = "Title of the plot"
    data_path = "path/to/plt-files"
    output_path = "path/to/save/results"
    variable_name = "mag_vorticity" #Variable to plot 

    # Run plotting
    plot_variable_slices(data_path, output_path, variable_name, case_title)

    # Create GIF
    gif_save_path = os.path.join(output_path, f"{variable_name}.gif")
    create_gif_from_images(output_path, gif_save_path)

```
