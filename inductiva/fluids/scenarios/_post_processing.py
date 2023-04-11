"""Post process SPH simulation outputs."""
import os

from base64 import b64encode
from IPython.display import HTML

from inductiva.fluids.post_processing import render_vtk
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
               fps: int = 10,
               color: str = "blue"):
        """Render the simulation as a movie.

        Args:
        fps: Number of frames per second to use in the movie. Renders a
            subset of the vtk files to create the movie.
            Default: 10.
        box_domain: Description of domain boundaries. As of now, it is fixed
            for a box defined by [x_min, x_max, y_min, y_max, z_min, z_max],
            but we can iterate to include objects ending in .obj, .stl or .vtk.
        color: The color of the markers in the simulation.
          If None, the default color is used.
        """

        vtk_dir = os.path.join(self.sim_output_dir, "vtk")

        movie_path = os.path.join(self.sim_output_dir, "movie.mp4")

        render_vtk(self.sim_output_dir,
                   movie_path,
                   camera=[(2., 2., 1.5), (0., 0., 0.), (1., 1., 3.)],
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
