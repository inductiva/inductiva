"""Post-processing and visualization utilities of the heat sink scenario."""

import os
import pathlib
from typing import Sequence

import pyvista as pv

from inductiva.types import Path
from inductiva.utils import files


class HeatSinkOutput:
    """Heat sink simulation output."""

    def __init__(self, sim_output_path: Path):
        """Initializes a `HeatSinkOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
        """

        self.sim_output_path = sim_output_path

    def render(
            self,
            movie_path: Path = "movie.mp4",
            fps: int = 10,
            cmap: str = "viridis",
            clim: Sequence[float] = (280, 300),
    ):
        """Renders temporal evolution of the temperature.
        
        A movie is produced representing the temporal evolution of the
        temperature in the heat sink and in the surrounding air flow.
         
        The movie is saved in the `movie_path` location.

        Args:
            movie_path: Path to the movie file.
            fps: Number of frames per second to use in the movie.
            cmap: Colormap used to represent temperature.
            clim: Colorbar limits.
        """

        movie_path = files.resolve_path(movie_path)

        # The OpenFOAM data reader from PyVista requires that a file named
        # "foam.foam" exists in the simulation output directory.
        # Create this file if it does not exist.
        foam_file_path = os.path.join(self.sim_output_path, "foam.foam")
        pathlib.Path(foam_file_path).touch(exist_ok=True)

        reader = pv.OpenFOAMReader(foam_file_path)

        plotter = pv.Plotter(off_screen=True)

        # Set camera position for a nice view.
        plotter.view_vector([-0.67, 0.49, -0.58], viewup=[0.31, 0.90, 0.29])
        plotter.camera.zoom(1.4)

        plotter.open_movie(movie_path, framerate=fps)

        for idx in range(1, reader.number_time_points):
            reader.set_active_time_point(idx)

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

            # Create custom color bar for temperature.
            plotter.add_scalar_bar(
                "Temperature [K]",
                position_x=0.45,
                width=0.5,
                title_font_size=18,
                label_font_size=18,
                fmt="%.0f",
                font_family="arial",
            )

            # Add time label.
            plotter.add_title(f"Time = {reader.active_time_value} [s]",
                              font_size=12)

            plotter.show_axes()

            plotter.write_frame()

        plotter.close()
