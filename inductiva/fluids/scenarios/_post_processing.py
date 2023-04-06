"""Post process SPH simulation outputs."""
import os
from typing import List, Optional

from base64 import b64encode
from IPython.display import HTML

import pyvista as pv

from inductiva.fluids.post_processing import render_vtk
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
               movie_fps: int = 10,
               box_domain: Optional[List[float]] = [0., 1., 0., 1., 0., 1.],
               color: str = "#00CCFF"):
        """Render the simulation as a movie.

        Args:
        movie_fps: The frames per second (fps) to be used in the movie.
          Default is 10 fps.
        box_domain: Description of domain boundaries. As of now, it is fixed
            for a box defined by [x_min, x_max, y_min, y_max, z_min, z_max],
            but we can iterate to include objects ending in .obj, .stl or .vtk.
        color: The color of the markers in the simulation.
          If None, the default color is used.
        """

        domain = pv.Box(bounds=box_domain)
        movie_path = os.path.join(self.sim_output_dir, "movie.mp4")

        render_vtk(self.sim_output_dir,
                                 movie_path,
                                 objects=domain,
                                 color=color,
                                 fps=movie_fps)

        with open(movie_path, "rb") as file_path:
            mp4 = file_path.read()
        movie_url = "data:video/mp4;base64," + b64encode(mp4).decode()

        return HTML(f"""
            <video alt="test" controls>
                <source src="{movie_url}" type="video/mp4">
            </video>
        """)
