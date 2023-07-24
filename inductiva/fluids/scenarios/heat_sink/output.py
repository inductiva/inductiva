"""TODO
"""
import os
from IPython.display import HTML

import pyvista as pv

from inductiva.types import Path
from inductiva.utils import files


class HeatSinkOutput:
    """Post-Process WindTunnel simulation outputs.

    Current Support:
        OpenFOAM
    """

    def __init__(self, sim_output_path: Path):
        """Initializes a `WindTunnelSimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
            time_step: Time step where we read the data.
        """

        self.sim_output_path = sim_output_path

    def render(
        self,
        movie_path: Path = "movie.mp4",
        fps: int = 10,
        cmap: str = "viridis",
        clim: list = [280, 300],
    ):
        """Render flow property over the object in the WindTunnel."""

        movie_path = files.resolve_path(movie_path)

        foam_file_path = os.path.join(self.sim_output_path, "foam.foam")

        # Create reading file
        with open(foam_file_path, "w", encoding="utf-8"):
            reader = pv.OpenFOAMReader(foam_file_path)

        plotter = pv.Plotter(off_screen=True)

        # Set camera position to a nice view.
        plotter.view_vector([-0.67, 0.49, -0.58], viewup=[0.31, 0.90, 0.29])
        plotter.camera.zoom(1.4)

        plotter.open_movie(movie_path, framerate=fps)

        for i in range(1, reader.number_time_points):
            reader.set_active_time_point(i)

            mesh = reader.read()

            fins_data = mesh["fins"]["boundary"]
            fluid_data = mesh["fluid"]["internalMesh"]

            fluid_data_slice = fluid_data.slice(normal=(-1, 0, 0))

            plotter.add_mesh(fins_data,
                             scalars="T",
                             cmap=cmap,
                             clim=clim,
                             show_scalar_bar=False)

            plotter.add_mesh(fluid_data_slice,
                             scalars="T",
                             cmap=cmap,
                             clim=clim,
                             show_scalar_bar=False)

            plotter.add_scalar_bar(
                "Temperature [K]",
                position_x=0.45,
                width=0.5,
                title_font_size=18,
                label_font_size=18,
                n_labels=5,
                fmt="%.0f",
                font_family="arial",
            )

            plotter.add_title(f"Time = {reader.active_time_value} [s]",
                              font_size=12)

            plotter.show_axes()

            plotter.write_frame()

        plotter.close()
