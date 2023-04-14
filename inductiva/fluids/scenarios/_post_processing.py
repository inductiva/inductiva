"""Post process SPH simulation outputs."""
import os

from base64 import b64encode
from IPython.display import HTML

from inductiva.utils.visualization import create_movie_from_vtk
from inductiva.types import Path


class SPHSimulationOutput:
    """Post process SPH simulation outputs."""

    def __init__(self, sim_output_path: Path):
        """Initializes a `SimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
            """
        self.sim_output_dir = sim_output_path

    def render(self,
               virtual_display: bool = True,
               fps: int = 10,
               color: str = "blue"):
        """Render the simulation as a movie.

        Args:
        fps: Number of frames per second to use in the movie. Renders a
            subset of the vtk files to create the movie.
            Default: 10.
        color: The color of the markers in the simulation.
          If None, the default color is used.
        """

        vtk_dir = os.path.join(self.sim_output_dir, "vtk")

        movie_path = os.path.join(self.sim_output_dir, "movie.mp4")

        create_movie_from_vtk(vtk_dir,
                              movie_path,
                              virtual_display=virtual_display,
                              camera=[(3., 3., 2.), (0., 0., 0.), (1., 1., 2.)],
                              fps=fps,
                              color=color)

        with open(movie_path, "rb") as file_path:
            mp4 = file_path.read()
        movie_url = "data:video/mp4;base64," + b64encode(mp4).decode()

        return HTML(f"""
            <video alt="test" controls>
                <source src="{movie_url}" type="video/mp4">
            </video>
        """)
