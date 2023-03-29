"""Visualization utilities."""

import os
import tempfile
from typing import Dict, List, Optional

from absl import logging

import imageio
import matplotlib
import matplotlib.colors as clr
import matplotlib.pyplot as plt
from tqdm import tqdm
import xarray as xr

MPL_CONFIG_PARAMS = {
    "font.size": 14,
    "axes.titlesize": "medium",
    "figure.dpi": 100,
}


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


def create_2d_scatter_plot(
    dataset: xr.Dataset,
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

        dataset.plot.scatter(
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


def create_3d_scatter_plot(
    dataset: xr.Dataset,
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

        dataset.plot.scatter(
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


def create_2d_scatter_plot_movie(
    dataset: xr.Dataset,
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
        var_iterator = dataset.groupby(iter_var)
        for i, (_, dataset_i) in tqdm(enumerate(var_iterator),
                                      total=len(var_iterator)):
            create_2d_scatter_plot(
                dataset=dataset_i,
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


def create_3d_scatter_plot_movie(
    dataset: xr.Dataset,
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
        var_iterator = dataset.groupby(iter_var)
        for i, (_, dataset_i) in tqdm(enumerate(var_iterator),
                                      total=len(var_iterator)):
            create_3d_scatter_plot(
                dataset=dataset_i,
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


def create_color_plot(
    data_array: xr.DataArray,
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

        data_array.plot.pcolormesh(
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


def create_color_plot_movie(
    data_array: xr.DataArray,
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
        var_iterator = data_array.groupby(iter_var)
        for i, (_, data_array_i) in tqdm(enumerate(var_iterator),
                                         total=len(var_iterator)):
            create_color_plot(
                data_array=data_array_i,
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
