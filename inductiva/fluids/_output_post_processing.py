"""Post process SPH simulation outputs."""
import os

from base64 import b64encode
from IPython.display import HTML

import xarray as xr

from inductiva.utils import visualization
from inductiva.types import Path


class SimulationOutput:
    """Post process SPH simulation outputs."""

    def __init__(self, sim_output_path: Path):
        """Initializes a `SimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
            """
        self.sim_output_dir = sim_output_path

    def render(self,
               movie_fps: int = 60,
               color: str = None,
               marker_size: float = 20.0,
               alpha: float = 1.0):
        """Render the simulation as a movie.

        Args:
        movie_fps: The frames per second (fps) to be used in the movie.
          Default is 60 fps.
        color: The color of the markers in the simulation.
          If None, the default color is used.
        marker_size: The size of the markers in the simulation. Default is 20.0.
        alpha: The opacity of the markers in the simulation.
          Default is 1.0 (fully opaque)
        """

        # Read simulation particle data
        particle_data = xr.open_mfdataset(
            os.path.join(self.sim_output_dir, "netcdf", "*.nc"))

        movie_path = os.path.join(self.sim_output_dir, "movie.mp4")

        visualization.create_3d_scatter_plot_movie(dataset=particle_data,
                                                   iter_var="time",
                                                   x_var="x",
                                                   y_var="y",
                                                   z_var="z",
                                                   movie_path=movie_path,
                                                   movie_fps=movie_fps,
                                                   x_limits=[0., 1.],
                                                   y_limits=[0., 1.],
                                                   z_limits=[0., 1.],
                                                   color=color,
                                                   marker_size=marker_size,
                                                   alpha=alpha)

        with open(movie_path, "rb") as file_path:
            mp4 = file_path.read()
        movie_url = "data:video/mp4;base64," + b64encode(mp4).decode()

        return HTML(f"""
            <video alt="test" controls>
                <source src="{movie_url}" type="video/mp4">
            </video>
        """)
