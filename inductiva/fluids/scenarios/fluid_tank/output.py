"""Post-processing and visualization utilities of the fluid tank scenario."""

import os
from typing import Optional

try:
    import xarray as xr
except ImportError:
    xr = None

from inductiva.types import Path
from inductiva.fluids.post_processing.splishsplash import convert_vtk_data_dir_to_netcdf
from inductiva.utils import files, visualization, optional_deps


class FluidTankOutput:
    """Fluid tank simulation output."""

    def __init__(self, sim_output_path: Path, output_time_step: float = 1):
        """Initializes a `FluidTankOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
        """

        self.sim_output_path = sim_output_path
        self.output_time_step = output_time_step

    @optional_deps.needs_fluids_extra_deps
    def render(
        self,
        movie_path: Path = "movie.mp4",
        fps: int = 10,
        color: Optional[str] = None,
    ):
        """Renders temporal evolution of the fluid particle positions.
        
        A movie is produced representing the temporal evolution of the positions
        of the particles representing the fluid in the tank. Each frame of the
        movie is a scatter plot of the particle positions at a given time.
         
        The movie is saved in the `movie_path` location.

        Args:
            movie_path: Path to the movie file.
            fps: Number of frames per second to use in the movie.
            color: Color used to represent the particles.
        """

        convert_vtk_data_dir_to_netcdf(
            data_dir=os.path.join(self.sim_output_path, "vtk"),
            output_time_step=self.output_time_step,
            netcdf_data_dir=os.path.join(self.sim_output_path, "netcdf"))

        particle_data = xr.open_mfdataset(
            os.path.join(self.sim_output_path, "netcdf", "*.nc"))

        movie_path = files.resolve_path(movie_path)

        visualization.create_3d_scatter_plot_movie(
            particle_data,
            iter_var="time",
            x_var="x",
            y_var="y",
            z_var="z",
            x_limits=[-0.5, 0.5],
            y_limits=[-0.5, 0.5],
            z_limits=[-0.1, 1],
            color=color,
            movie_path=movie_path,
            movie_fps=fps,
        )
