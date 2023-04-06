"""Visualization for fluid dynamics from .vtk files."""

import os

import pyvista as pv

from inductiva.fluids.post_processing import process_vtk


def render_vtk(sim_output_dir: str,
               movie_path: str,
               objects: pv.PolyData = None,
               color: str = "#00CCFF",
               point_size: float = 5.0,
               point_opacity: float = 1.0,
               fps: int = 10) -> None:
    """Creates movie from a series of vtk files.

    The order of the vtk file name determines the order with which they
    are render in the movie. For example, vtk file 'particle_001.vtk' will
    appear before 'particle_002.vtk'.

    Args:
        sim_output_dir: Directory containing a 'vtk' folder.
        objects: Object of pyvista.PolyData type describing the domain or
            an object inside.
        color: The color of the points in the simulation.
          The default color is light blue.
        point_size: The size of the points in the simulation. Default is 5.0.
        point_opacity: The opacity of the markers in the simulation.
          Default is 1.0 (fully opaque)
        fps: Number of frames per second to use in the movie. This cuts some
            frames of data for speed and lower quality purposes.
    """
    vtk_dir = os.path.join(sim_output_dir, "vtk")

    vtk_files = process_vtk.get_sorted_vtk_files(vtk_dir)

    plotter = pv.Plotter()
    fps=30
    # Open a movie file
    plotter.open_movie(movie_path, framerate=fps)

    # Set camera
    camera = pv.Camera()
    plotter.camera.position = (2.5, 2.5, 1.5)
    plotter.camera.focal_point=(0.,0.,0.)

    # Set domain surroundings
    if objects is not None:
        plotter.add_mesh(objects, opacity=0.15)

    # Set initial frame that defines the following ones
    initial_vtk_file = os.path.join(sim_output_dir, vtk_files[0])
    mesh = pv.read(initial_vtk_file)
    plotter.add_mesh(mesh,
                     color=color,
                     point_size=point_size,
                     opacity=point_opacity,
                     render_points_as_spheres=True)

    plotter.write_frame()

    for index, frame_file in enumerate(vtk_files):
        if index % int(round(60/fps)) == 0:
            frame_path = os.path.join(sim_output_dir, frame_file)
            frame = pv.read(frame_path)
            mesh.copy_from(frame)
            plotter.render()
            plotter.write_frame()

    plotter.close()
