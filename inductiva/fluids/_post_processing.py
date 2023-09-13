"""Post process SPH simulation outputs."""
import os

from inductiva.utils import visualization, optional_deps
from inductiva.types import Path


class SPHSimulationOutput:
    """Post process SPH simulation outputs."""

    def __init__(self, sim_output_path: Path):
        """Initializes a `SimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
            """
        self.sim_output_dir = sim_output_path

    @optional_deps.needs_fluids_extra_deps
    def render(self,
               virtual_display: bool = False,
               movie_path: str = "movie.mp4",
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

        visualization.create_movie_from_vtk(vtk_dir,
                                            movie_path,
                                            virtual_display=virtual_display,
                                            camera=[(3., 3., 2.), (0., 0., 0.),
                                                    (1., 1., 2.)],
                                            fps=fps,
                                            color=color)
