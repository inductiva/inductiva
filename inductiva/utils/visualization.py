"""Visualization utilities."""
# pylint: disable=protected-access

import os
import tempfile
from typing import Dict, List, Optional
import base64
import io
from time import sleep

from absl import logging
from tqdm import tqdm
try:
    import imageio
    import matplotlib
    import matplotlib.colors as clr
    import matplotlib.pyplot as plt
except ImportError:
    imageio = None
    matplotlib = None
    clr = None
    plt = None

try:
    import ipywidgets
    import PIL
except ImportError:
    ipywidgets = None
    PIL = None

try:
    import pyvista as pv
    import xarray as xr
except ImportError:
    pv = None
    xr = None

from inductiva.utils import files, optional_deps
import threading

MPL_CONFIG_PARAMS = {
    "font.size": 14,
    "axes.titlesize": "medium",
    "figure.dpi": 100,
}


@optional_deps.needs_molecules_extra_deps
def create_movie_from_widget(view,
                             output_path="movie.mp4",
                             fps=8,
                             start=0,
                             stop=-1,
                             step=1,
                             timeout=0.2):
    """Create a movie from a ipywidget view.
    Args: 
        view: ipywidget visualization.
        output_path: name of the output movie file.
        fps: frame rate (frames per second). 
        start: starting frame. 
        stop: ending frame number.
        step: step between frames.
        timeout: the waiting time between rendering two consecutive frames.
    """
    if stop < 0:
        stop = view.max_frame + 1
    frames_range = range(start, stop, step)
    max_frame_digits = len(str(stop))
    event = threading.Event()

    def _make(event):
        with tempfile.TemporaryDirectory() as tmp_dir:
            for i in tqdm(frames_range):
                if not event.is_set():
                    view.frame = i
                    sleep(timeout)  # time to update the view
                    iw = view.render_image()
                    image_data = view._image_data
                    sleep(timeout)
                    filename = os.path.join(
                        tmp_dir, f"frame-{str(i).zfill(max_frame_digits)}.png")
                    try:
                        decode_save_image(image_data, filename)
                    except PIL.UnidentifiedImageError:
                        print(f"Error: Unidentified image at frame {i}")
                        continue
                    iw.close()
                    sleep(timeout)

            if not event.is_set():
                with ipywidgets.Output():
                    create_movie_from_frames(tmp_dir, output_path, fps)

    thread = threading.Thread(target=_make, args=(event,))
    thread.daemon = True
    thread.start()


def decode_save_image(image_data, filename):
    """Decode image data bytes and save to file."""
    im_bytes = base64.b64decode(image_data)
    im_bytes = io.BytesIO(im_bytes)
    image = PIL.Image.open(im_bytes)
    image.save(filename, "PNG")


def create_movie_from_frames(frames_dir: str,
                             movie_path: str,
                             fps: int = 10) -> None:
    """Creates movie from a series of png image files.

    The order of the png image file names determines the order with which they
    are rendered in the movie. For example, image 'frame-001.png' will appear
    before 'frame-002.png'.

    Args:
        frames_dir: Directory where the png files are stored.
        movie_path: Path to the movie to be created.
        fps: Number of frames per second to use in the movie.
    """

    writer = imageio.get_writer(movie_path, fps=fps)

    frames_list = [
        file for file in os.listdir(frames_dir) if file.endswith(".png")
    ]
    frames_list.sort()

    for file_name in frames_list:
        frame_path = os.path.join(frames_dir, file_name)
        writer.append_data(imageio.imread(frame_path))

    writer.close()


def create_movie_from_vtk(vtk_output_dir: str,
                          movie_path: str,
                          virtual_display: bool = True,
                          scalars: str = None,
                          scalar_limits: Optional[List[float]] = None,
                          camera=None,
                          color: str = "blue",
                          cmap: str = None,
                          fps: int = 10) -> None:
    """Creates movie from a series of vtk files.

    The order of the vtk file name determines the order with which they
    are rendered in the movie. For example, vtk file 'frame_001.vtk' will
    appear before 'frame_002.vtk'.

    Args:
        vtk_output_dir: Directory containing the vtk files.
        movie_path: Path to save the movie.
        virtual_display: Whether to use a virtual display to render
            the movie.
        scalar: Scalars used to “color” the mesh. Accepts a string name
            of an array that is present on the mesh or an array equal to
            the number of cells or the number of points in the mesh.
            Array should be sized as a single vector.
        scalar_bounds: Color bar range for scalars. Defaults to minimum
            and maximum of scalars array. Example: [-1, 2].
        objects: Object of pyvista.PolyData type describing the domain or
            an object inside.
        camera: Camera description must be one of the following:
          - List of three tuples describing the position, focal-point
          and view-up: [(2.0, 5.0, 13.0), (0.0, 0.0, 0.0), (-0.7, -0.5, 0.3)]
          - List with a view-vector: [-1.0, 2.0, -5.0]
          - A string with the plane orthogonal to the view direction: 'xy'.
          https://docs.pyvista.org/api/plotting/_autosummary/pyvista.CameraPosition.html
        color: Color of the points used to represent particles. Default: "blue".
        cmap: string with the name of the matplotlib colormap to use
            when mapping the scalars. See available Matplotlib colormaps.
        fps: Number of frames per second to use in the movie. Renders only a
            subset of the vtk files to create the movie. This is done for
            speed purposes.
            Default: 10.
    """
    if virtual_display:
        pv.start_xvfb()

    pv.global_theme.background = "white"
    vtk_output_dir = files.resolve_path(vtk_output_dir)
    vtk_files = files.get_sorted_files(vtk_output_dir, ".vtk")

    with tempfile.TemporaryDirectory() as tmp_dir:
        logging.info("Creating movie frames...")
        for index, frame_file in enumerate(vtk_files):
            if index % int(round(60 / fps)) == 0:
                frame_path = os.path.join(vtk_output_dir, frame_file)
                image_frame_path = os.path.join(
                    tmp_dir, "frame_" + str(index).zfill(5) + ".png")

                create_frame_from_vtk(frame_path,
                                      image_frame_path=image_frame_path,
                                      camera=camera,
                                      scalars=scalars,
                                      scalar_limits=scalar_limits,
                                      color=color,
                                      cmap=cmap)

        logging.info("Creating movie '%s'.", movie_path)
        create_movie_from_frames(frames_dir=tmp_dir,
                                 movie_path=movie_path,
                                 fps=fps)


def create_frame_from_vtk(frame_path: str,
                          image_frame_path: str,
                          scalars: str = None,
                          scalar_limits: Optional[List[float]] = None,
                          camera=None,
                          color: str = None,
                          cmap: str = None):
    """Render a .png image from a vtk file."""

    frame = pv.read(frame_path)

    frame.plot(off_screen=True,
               screenshot=image_frame_path,
               cpos=camera,
               scalars=scalars,
               clim=scalar_limits,
               render_points_as_spheres=True,
               color=color,
               cmap=cmap)


@optional_deps.needs_fluids_extra_deps
def create_2d_scatter_plot(
    xr_dataset: "xr.Dataset",
    x_var: str,
    y_var: str,
    image_path: str,
    marker_size: float = 20,
    alpha: float = 1,
    color: Optional[str] = None,
    color_var: Optional[str] = None,
    x_limits: Optional[List[float]] = None,
    y_limits: Optional[List[float]] = None,
    color_limits: Optional[List[float]] = None,
    color_map: Optional[str] = "viridis",
    aspect: Optional[str] = None,
    matplotlib_config_params: Optional[Dict] = None,
):
    """Creates a 2d scatter plot of a dataset."""

    if matplotlib_config_params is None:
        matplotlib_config_params = MPL_CONFIG_PARAMS

    with matplotlib.rc_context(rc=matplotlib_config_params):
        fig = plt.figure()
        ax = fig.add_subplot()

        color_norm = None
        if color_var is not None and color_limits is not None:
            color_norm = clr.Normalize(*color_limits)

        xr_dataset.plot.scatter(
            ax=ax,
            x=x_var,
            y=y_var,
            hue=color_var,
            s=marker_size,
            alpha=alpha,
            c=color,
            xlim=x_limits,
            ylim=y_limits,
            norm=color_norm,
            cmap=color_map,
        )

        if aspect is not None:
            ax.set(aspect=aspect)

        fig.tight_layout()
        fig.savefig(image_path)
        plt.close(fig)


@optional_deps.needs_fluids_extra_deps
def create_3d_scatter_plot(
    xr_dataset: "xr.Dataset",
    x_var: str,
    y_var: str,
    z_var: str,
    image_path: str,
    marker_size: float = 20,
    alpha: float = 1,
    color: Optional[str] = None,
    color_var: Optional[str] = None,
    x_limits: Optional[List[float]] = None,
    y_limits: Optional[List[float]] = None,
    z_limits: Optional[List[float]] = None,
    color_limits: Optional[List[float]] = None,
    color_map: Optional[str] = "viridis",
    matplotlib_config_params: Optional[Dict] = None,
):
    """Creates a 3d scatter plot of a dataset."""

    if matplotlib_config_params is None:
        matplotlib_config_params = MPL_CONFIG_PARAMS

    with matplotlib.rc_context(rc=matplotlib_config_params):
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")

        color_norm = None
        if color_var is not None and color_limits is not None:
            color_norm = clr.Normalize(*color_limits)

        xr_dataset.plot.scatter(
            ax=ax,
            x=x_var,
            y=z_var,
            z=y_var,
            s=marker_size,
            alpha=alpha,
            c=color,
            hue=color_var,
            xlim=x_limits,
            ylim=z_limits,
            norm=color_norm,
            cmap=color_map,
        )

        ax.set_zlim(y_limits)

        fig.tight_layout()
        fig.savefig(image_path)
        plt.close(fig)


@optional_deps.needs_fluids_extra_deps
def create_2d_scatter_plot_movie(
    xr_dataset: "xr.Dataset",
    iter_var: str,
    x_var: str,
    y_var: str,
    movie_path: str,
    movie_fps: int = 10,
    marker_size: float = 20,
    alpha: float = 1,
    color: Optional[str] = None,
    color_var: Optional[str] = None,
    x_limits: Optional[List[float]] = None,
    y_limits: Optional[List[float]] = None,
    color_limits: Optional[List[float]] = None,
    color_map: Optional[str] = "viridis",
    aspect: Optional[str] = None,
    matplotlib_config_params: Optional[Dict] = None,
):
    """Creates a movie of the dataset with 2d scatter plots."""

    with tempfile.TemporaryDirectory() as tmp_dir:
        logging.info("Creating movie frames...")
        var_iterator = xr_dataset.groupby(iter_var)
        for i, (_, dataset_i) in tqdm(enumerate(var_iterator),
                                      total=len(var_iterator)):
            create_2d_scatter_plot(
                xr_dataset=dataset_i,
                x_var=x_var,
                y_var=y_var,
                marker_size=marker_size,
                alpha=alpha,
                color=color,
                color_var=color_var,
                x_limits=x_limits,
                y_limits=y_limits,
                color_limits=color_limits,
                color_map=color_map,
                aspect=aspect,
                image_path=os.path.join(tmp_dir, f"frame_{i:06d}.png"),
                matplotlib_config_params=matplotlib_config_params,
            )

        logging.info("Finished creating movie frames.")

        logging.info("Creating movie '%s'.", movie_path)
        create_movie_from_frames(frames_dir=tmp_dir,
                                 movie_path=movie_path,
                                 fps=movie_fps)


@optional_deps.needs_fluids_extra_deps
def create_3d_scatter_plot_movie(
    xr_dataset: "xr.Dataset",
    iter_var: str,
    x_var: str,
    y_var: str,
    z_var: str,
    movie_path: str,
    movie_fps: int = 10,
    marker_size: float = 20,
    alpha: float = 1,
    color: Optional[str] = None,
    color_var: Optional[str] = None,
    x_limits: Optional[List[float]] = None,
    y_limits: Optional[List[float]] = None,
    z_limits: Optional[List[float]] = None,
    color_limits: Optional[List[float]] = None,
    color_map: Optional[str] = "viridis",
    matplotlib_config_params: Optional[Dict] = None,
):
    """Creates a movie of the dataset with 3d scatter plots."""

    with tempfile.TemporaryDirectory() as tmp_dir:
        logging.info("Creating movie frames...")
        var_iterator = xr_dataset.groupby(iter_var)
        for i, (_, dataset_i) in tqdm(enumerate(var_iterator),
                                      total=len(var_iterator)):
            create_3d_scatter_plot(
                xr_dataset=dataset_i,
                x_var=x_var,
                y_var=y_var,
                z_var=z_var,
                marker_size=marker_size,
                alpha=alpha,
                color=color,
                color_var=color_var,
                x_limits=x_limits,
                y_limits=y_limits,
                z_limits=z_limits,
                color_limits=color_limits,
                color_map=color_map,
                image_path=os.path.join(tmp_dir, f"frame_{i:06d}.png"),
                matplotlib_config_params=matplotlib_config_params,
            )

        logging.info("Finished creating movie frames.")

        logging.info("Creating movie '%s'.", movie_path)
        create_movie_from_frames(frames_dir=tmp_dir,
                                 movie_path=movie_path,
                                 fps=movie_fps)


@optional_deps.needs_fluids_extra_deps
def create_color_plot(
    xr_data_array: "xr.DataArray",
    image_path: str,
    x_limits: Optional[List[float]] = None,
    y_limits: Optional[List[float]] = None,
    color_limits: Optional[List[float]] = None,
    color_map: Optional[str] = "viridis",
    aspect: Optional[str] = None,
    matplotlib_config_params: Optional[Dict] = None,
):
    """Creates a color plot of a data array."""

    if matplotlib_config_params is None:
        matplotlib_config_params = MPL_CONFIG_PARAMS

    with matplotlib.rc_context(rc=matplotlib_config_params):
        fig = plt.figure()
        ax = fig.add_subplot()

        color_norm = None
        if color_limits is not None:
            color_norm = clr.Normalize(*color_limits)

        xr_data_array.plot.pcolormesh(
            ax=ax,
            xlim=x_limits,
            ylim=y_limits,
            norm=color_norm,
            cmap=color_map,
        )

        if aspect is not None:
            ax.set(aspect=aspect)

        fig.tight_layout()
        fig.savefig(image_path)
        plt.close(fig)


@optional_deps.needs_fluids_extra_deps
def create_color_plot_movie(
    xr_data_array: "xr.DataArray",
    iter_var: str,
    movie_path: str,
    movie_fps: int = 10,
    x_limits: Optional[List[float]] = None,
    y_limits: Optional[List[float]] = None,
    color_limits: Optional[List[float]] = None,
    color_map: Optional[str] = "viridis",
    aspect: Optional[str] = None,
    matplotlib_config_params: Optional[Dict] = None,
):
    """Creates a movie of the data array with color plots."""

    with tempfile.TemporaryDirectory() as tmp_dir:
        logging.info("Creating movie frames...")
        var_iterator = xr_data_array.groupby(iter_var)
        for i, (_, data_array_i) in tqdm(enumerate(var_iterator),
                                         total=len(var_iterator)):
            create_color_plot(
                xr_data_array=data_array_i,
                x_limits=x_limits,
                y_limits=y_limits,
                color_limits=color_limits,
                color_map=color_map,
                aspect=aspect,
                image_path=os.path.join(tmp_dir, f"frame_{i:06d}.png"),
                matplotlib_config_params=matplotlib_config_params,
            )

        logging.info("Finished creating movie frames.")

        logging.info("Creating movie '%s'.", movie_path)
        create_movie_from_frames(frames_dir=tmp_dir,
                                 movie_path=movie_path,
                                 fps=movie_fps)
