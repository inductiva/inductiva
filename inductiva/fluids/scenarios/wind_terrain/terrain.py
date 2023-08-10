"""Represent a complex terrain profile."""
import typing
import random

import pyvista as pv

from inductiva.generative import procedural


class Terrain:
    """Represents a Terrain profile.
    
    A Terrain profile is defined by the surface of the terrain.
    Here, we describe it by a mesh with coordinates (x, y, z).
    The (x, y) reflect the (longitude, latitude) of the terrain
    and the z coordinate is the elevation of the terrain.
    """

    def __init__(self, terrain_mesh):
        """Initializes the `Terrain` object.
        
        Args:
            terrain_mesh (pyvista.StructuredGrid): Terrain mesh."""
        self.mesh = terrain_mesh

    @classmethod
    def from_file(cls, terrain_file: str):
        """Creates a `Terrain` object from a text file."""
        terrain_mesh = pv.read(terrain_file)

        return cls(terrain_mesh)

    @classmethod
    def from_random_generation(cls,
                               x_range: typing.Sequence[float],
                               y_range: typing.Sequence[float],
                               x_num: int,
                               y_num: int,
                               height_factor: float = 10,
                               initial_roughness: float = 1,
                               roughness_factor: float = 0.5):
        """Creates a `Terrain` object with random elevations.
        The elevation of the corners are randomly selected within
        the range of 0 to 1. This limits the maximum terrain height
        due to the nature of the terrain generation algorithm.
        To allow for more complex terrains, we allow the user to multiply
        by a `heigh_factor` to increase the maximum height of the terrain.
        In this way, the maximum_height possible will be around 
        the value of the `height_factor`.

        TODO: Improve terrain generation for high discrepancy
        terrains.

        Args:
            x_range = [x_min, x_max]: Range of the x-axis in meters (m).
            y_range = [y_min, y_max]: Range of the y-axis in meters (m).
            x_num: Number of points in the x-axis.
            y_num: Number of points in the y-axis.
            height_factor: Factor to multiply the elevation by.
            initial_roughness: Initial roughness of the terrain.
            roughness_factor: Factor to multiply the roughness by.
        """
        corner_values = [
            random.uniform(0, 1),
            random.uniform(0, 1),
            random.uniform(0, 1),
            random.uniform(0, 1)
        ]

        x_grid, y_grid, z_elevation = procedural.generate_random_terrain(
            x_range, y_range, x_num, y_num, corner_values, initial_roughness,
            roughness_factor)

        terrain = pv.StructuredGrid(x_grid, y_grid, z_elevation * height_factor)

        return cls(terrain)

    def to_text_file(self, text_file_path: str):
        """Saves the terrain to a text file.
        
        Args
            text_file_path: Path to the text file
        """

        terrain_geometry = self.mesh.extract_geometry()
        terrain_geometry.save(text_file_path)

    @property
    def bounds(self):
        """Returns the bounds of the terrain."""

        terrain_bounds = {
            "x": [self.mesh.bounds[0], self.mesh.bounds[1]],
            "y": [self.mesh.bounds[2], self.mesh.bounds[3]],
            "z": [self.mesh.bounds[4], self.mesh.bounds[5]]
        }

        return terrain_bounds

    @property
    def center(self):
        """Returns the center of the terrain."""

        terrain_center = {
            "x": (self.mesh.bounds[0] + self.mesh.bounds[1]) / 2,
            "y": (self.mesh.bounds[2] + self.mesh.bounds[3]) / 2
        }

        return terrain_center

    @property
    def mesh_resolution(self):
        """Returns the resolution of the terrain mesh."""

        num_cells = self.mesh.n_cells
        num_points = self.mesh.n_points

        return num_cells, num_points

    def plot(self,
             off_screen: bool = False,
             save_path: str = None,
             background_color: str = "white",
             colormap: str = "cividis",
             lighting: bool = True):
        """Renders the terrain.
        
        Args:
            off_screen: If True, the plot is not shown.
            save_path: Path to save the screenshot.
            background_color: Background color of the plot.
            colormap: Colormap to use.
            lighting: Enable or disable view direction lighting.
                It highlights the details in the terrain.
        TODO: camera play
        """

        # Create Elevation from z-values to plot by it.
        self.mesh["Elevation"] = self.mesh.z.T.reshape(-1)
        plotter = pv.Plotter(off_screen=off_screen)
        plotter.background_color = background_color
        plotter.add_mesh(self.mesh,
                         scalars="Elevation",
                         lighting=lighting,
                         cmap=colormap)
        plotter.show(screenshot=save_path)
        plotter.close()
