"""Visualization for fluid dynamics from .vtk files."""

import os
import tempfile
from typing import List, Optional

from absl import logging
import pyvista as pv

from inductiva.utils.visualization import create_movie_from_frames
from inductiva.fluids.post_processing import process_vtk

def render_vtk(vtk_output_dir: str,
               movie_path: str,
               scalars: str = None,
               scalar_bounds: Optional[List[float]] = None,
               camera=None,
               color: str = "blue",
               cmap: str = None,
               fps: int = 10) -> None:
    """Creates movie from a series of vtk files.

    The order of the vtk file name determines the order with which they
    are rendered in the movie. For example, vtk file 'particle_001.vtk' will
    appear before 'particle_002.vtk'.

    Args:
        vtk_output_dir: Directory containing the vtk files.
        movie_path: Path to save the movie.
        scalar: Scalars used to “color” the mesh. Accepts a string name
            of an array that is present on the mesh or an array equal to
            the number of cells or the number of points in the mesh.
            Array should be sized as a single vector.
        scalar_bounds: Color bar range for scalars. Defaults to minimum
            and maximum of scalars array. Example: [-1, 2].
        objects: Object of pyvista.PolyData type describing the domain or
            an object inside.
        camera: Camera description must be one of the following:
          - List of three tuples describing the position, focal-point
          and view-up: [(2.0, 5.0, 13.0), (0.0, 0.0, 0.0), (-0.7, -0.5, 0.3)]
          - List with a view-vector: [-1.0, 2.0, -5.0]
          - A string with the plane orthogonal to the view direction: 'xy'.
          https://docs.pyvista.org/api/plotting/_autosummary/pyvista.CameraPosition.html
        color: Color of the points used to represent particles. Default: "blue".
        cmap: string with the name of the matplotlib colormap to use
            when mapping the scalars. See available Matplotlib colormaps.
        fps: Number of frames per second to use in the movie. Renders only a
            subset of the vtk files to create the movie. This is done for
            speed purposes.
            Default: 10.
    """

    vtk_files = process_vtk.get_sorted_vtk_files(vtk_output_dir)

    with tempfile.TemporaryDirectory() as tmp_dir:
        logging.info("Creating movie frames...")
        for index, frame_file in enumerate(vtk_files):
            if index % int(round(60 / fps)) == 0:
                frame_path = os.path.join(vtk_output_dir, frame_file)
                image_frame_path = os.path.join(
                    tmp_dir, "Frame_" + str(index).zfill(5) + ".png")

                render_vtk_frame(frame_path,
                                 image_frame_path=image_frame_path,
                                 camera=camera,
                                 scalars=scalars,
                                 scalar_bounds=scalar_bounds,
                                 color=color,
                                 cmap=cmap)

        logging.info("Creating movie '%s'.", movie_path)
        create_movie_from_frames(frames_dir=tmp_dir,
                                 movie_path=movie_path,
                                 fps=fps)


def render_vtk_frame(frame_path: str,
                     image_frame_path: str,
                     scalars: str = None,
                     scalar_bounds: Optional[List[float]] = None,
                     camera=None,
                     color: str = None,
                     cmap: str = None):
    """Render a .png image from a vtk file."""

    frame = pv.read(frame_path)

    frame.plot(off_screen=True,
               screenshot=image_frame_path,
               cpos=camera,
               scalars=scalars,
               clim=scalar_bounds,
               render_points_as_spheres=True,
               color=color,
               cmap=cmap)
