# AMR-Wind Post-Processing

The easiest way to get started with post-processing AMR-Wind results is by using [yt](https://yt-project.org/doc/visualizing/plots.html#slice-plots), 
a Python-based package designed for analyzing and visualizing volumetric simulation data. 
It supports AMReX plotfiles out of the box and is ideal for scripting and automation. 
Alternatively, [Paraview](https://www.paraview.org/) also natively supports AMReX plotfiles and provides an interactive, GUI-based approach for quickly exploring and visualizing simulation results.

## Using `yt` for Post-Processing

yt is ideal for scripting and automated analysis of AMR data. It allows quick generation of slices, projections, and derived fields with just a few lines of Python code.

### 1. Install yt (if not installed already)

```bash
pip install yt
```

### 2. Basic yt script to load and plot a slice from the volume data

```python
import yt

ds = yt.load("plt01000")  # Load the time step of choosing from the provided path
slc = yt.SlicePlot(ds, "z", ("boxlib", "velocityx")) # Create slice plot of x velocity - normal to z-axis
slc.set_cmap(("boxlib", "velocityx"), "plasma") # Choosing colormap
slc.set_log(("boxlib", "velocityx"), False) # To opt linear scale- log scale is default
slc.annotate_title("Velocity field slice")
slc.save("xvelocity_plot.png")

```

### 3. Loop over time steps to plot slices

The below script can be used to loop over all `plt00000`files (time steps) and extract a data slice to save it as a 2D plot.

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
    slc.save() # will save the plot files to the same folder as the plt files

```

### Sample script- Create an animation from plt files

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
    vmin, vmax = 0, 50  # Minimum and maximum values for color scale

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
    variable_name = "mag_vorticity" #Variable to be plotted (eg: velocityx, velocityy, mag_vorticity)

    # Run plotting
    plot_variable_slices(data_path, output_path, variable_name, case_title)

    # Create GIF
    gif_save_path = os.path.join(output_path, f"{variable_name}.gif")
    create_gif_from_images(output_path, gif_save_path)

```