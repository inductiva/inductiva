"""Post process SPlisHSPlasH simulation outputs."""
import os

from base64 import b64encode
from IPython.display import HTML

import xarray as xr

from inductiva_data import visualization
from inductiva.types import Path


class SimulationOutput:
    """Post process SPlisHSPlasH simulation outputs."""

    def __init__(self, sim_output_path: Path) -> None:
        """Initializes a `SimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
            """
        self.sim_output_dir = sim_output_path

        

    def render(self):
        """Generate a simulation movie.

        Args:
            color_quantity: Quantity to represent in the color scale of the
                scatter plot."""

        # Read simulation particle data
        particle_data_dir = os.path.join(self.sim_output_dir, "netcdf")

        particle_data = xr.open_mfdataset(os.path.join(particle_data_dir, "*.nc"))

        movie_path = os.path.join(self.sim_output_dir, "movie.mp4")

        visualization.create_3d_scatter_movie(
            dataset=particle_data,
            iter_var="time",
            x_var="x",
            y_var="y",
            z_var="z",
            movie_path=movie_path)

        with open(movie_path, "rb") as file_path:
            mp4 = file_path.read()
        movie_url = "data:video/mp4;base64," + b64encode(mp4).decode()

        return HTML(f"""
            <video alt="test" controls>
                <source src="{movie_url}" type="video/mp4">
            </video>
        """)
