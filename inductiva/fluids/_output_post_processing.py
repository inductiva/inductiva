"""Post process SPlisHSPlasH simulation outputs."""
import os

from IPython.display import HTML
from base64 import b64encode
import xarray as xr

from inductiva_data.visualization import create_3d_scatter_plot_movie

from inductiva.types import Path


class SimulationOutput:
    """Post process SPlisHSPlasH simulation outputs."""

    def __init__(self, sim_output_path: Path) -> None:
        """Initializes a `SimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
            """
        self.sim_output_dir = sim_output_path

    def render(self, color_quantity: str = None):
        """Generate a simulation movie.

        Args:
            color_quantity: Quantity to represent in the color scale of the
                scatter plot."""

        particle_data = xr.open_mfdataset(
            os.path.join(self.sim_output_dir, "netcdf", "*.nc"))

        movie_path = os.path.join(self.sim_output_dir, "movie.mp4")
        create_3d_scatter_plot_movie(
            particle_data,
            iter_var="time",
            x_var="x",
            y_var="y",
            z_var="z",
            color_var=color_quantity,
            movie_path=movie_path,
        )

        with open(movie_path, "rb") as fp:
            mp4 = fp.read()
        movie_url = "data:video/mp4;base64," + b64encode(mp4).decode()

        return HTML(f"""
            <video alt="test" controls>
                <source src="{movie_url}" type="video/mp4">
            </video>
        """)
