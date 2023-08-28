"""Bathymetry class and utilities."""

import random
from typing import Optional, Sequence, Tuple, Union

import matplotlib
import numpy as np
import scipy
import utm

import inductiva


class Bathymetry:
    """Represents a bathymetric profile.
    
    A bathymetric profile defines the depth of the sea bottom as a function of
    space, here described in Cartesian coordinates (x, y).
    
    Here, a bathymetry is represented with:
    - a 1D array of depths, in meters, measured at arbitrary points in space.
      Positive depths are below the water level.
    - two 1D arrays representing the x and y coordinates of the points where the
      depths are defined, in meters.
    """

    def __init__(
        self,
        depths: np.ndarray,
        x: np.ndarray,
        y: np.ndarray,
    ):
        """Initializes a `Bathymetry` object.
        
        Args:
            depths: A 1D array with the depths, in meters.
            x: A 1D array with the x coordinates of the points where depths are
              defined, in meters.
            y: Same as `x`, but for the y coordinates.
        """
        self.depths = depths
        self.x = x
        self.y = y

    @classmethod
    def from_bot_file(
        cls,
        bot_file_path: str,
        x_range: Sequence[float],
        y_range: Sequence[float],
    ):
        """Creates a `Bathymetry` object from a text file.
        
        The depth values are read from a text file. The text file must contain
        a 2D array with the depths, in meters. The first and second dimensions
        of the array in the text file (i.e. rows and columns) correspond to the
        x and y directions, respectively.

        Args:
            text_file_path: Path to the text file.
            x_range: The range of x values, in meters.
            y_range: The range of y values, in meters.
        """

        depths = np.loadtxt(bot_file_path)

        x, y = np.meshgrid(np.linspace(*x_range, depths.shape[0]),
                           np.linspace(*y_range, depths.shape[1]),
                           indexing="ij")

        return cls(depths.flatten(), x.flatten(), y.flatten())

    @classmethod
    def from_ascii_xyz_file(
        cls,
        ascii_xyz_file_path: str,
        remove_offset: bool = True,
    ):
        """Creates a `Bathymetry` object from an ASCII XYZ file.
        
        ASCII XYZ files store bathymetric data in a table where each line
        corresponds to a location (latitude, longitude pair). Several columns
        are available to characterize the depth at each location, namely:
        - longitude, in decimal degrees, e.g. 52.07334567;
        - latitude, in decimal degrees, e.g. 3.06033283;
        - the minum depth, in meters, e.g. 35.81;
        - the maximum depth, in meters, e.g. 35.81;
        - average depth (over a set of measurements), in meters, e.g. 35.81.
        - the standard deviation of the depth (over a set of measurements), in
          meters, e.g. 35.81;

        For more information, see
        https://emodnet.ec.europa.eu/sites/emodnet.ec.europa.eu/files/public/20171127_DTM_exchange_format_specification_v1.6.pdf.

        In this method, we:
        - load the latitude and longitude values from the file, and convert
          them to (x, y) UTM coordinates (see
          https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system
          for more information).
        - use average depth as the local depth value.

        Args:
            ascii_xyz_file_path: Path to the ASCII XYZ file.
            remove_offset: Whether to remove the offset of the x and y
              coordinates. If `True`, the minimum x and y values are subtracted
              from all x and y values, respectively.
        """

        # Load only the first, second and fifth columns of the file,
        # corresponding to longitude, latitude and average depth, respectively.
        bathymetry_data = np.loadtxt(ascii_xyz_file_path,
                                     usecols=(0, 1, 4),
                                     delimiter=";")

        x_lon = bathymetry_data[:, 0]
        y_lat = bathymetry_data[:, 1]
        depths = bathymetry_data[:, 2]

        x, y, _, _ = utm.from_latlon(y_lat, x_lon)

        if remove_offset:
            x = x - x.min()
            y = y - y.min()

        return cls(depths, x, y)

    @classmethod
    def from_random_depths(
        cls,
        x_range: Sequence[float],
        y_range: Sequence[float],
        x_num: int,
        y_num: int,
        max_depth: float = 10,
        initial_roughness: float = 1,
        roughness_factor: float = 0.5,
        percentile_above_water: float = 20,
    ):
        """Creates a `Bathymetry` object with random depths.

        The depths of the corners of the grid are chosen according to a maximum
        depth value `max_depth` and a percentile of the domain that must be
        above water `percentile_above_water`. The corners on the lower x
        boundary (East) are assumed to be below water (i.e. have 0 < depths <
        `max_depth`). The corners on the upper x boundary (West) are assumed to
        be above water (i.e. have - `max_depth` < depths < 0).
        
        Args:
            x_range: The range of x values, in meters.
            y_range: The range of y values, in meters.
            x_num: Number of grid points in the x direction.
            y_num: Number of grid points in the y direction.
            max_depth: Maximum depth value, in meters.
            initial_roughness: Initial roughness value, in meters. Controls the
              initial range of randomness of the Diamond-Square algorithm.
            roughness_factor: Roughness factor. Must be between 0 and 1.
              Controls the rate at which the range of randomness of the
              Diamond-Square algorithm decreases over iterations.
            percentile_above_water: Percentile of the depths that must be above
              water. Must be between 0 and 100.
        """

        corner_values = [
            random.uniform(0, max_depth),
            random.uniform(0, max_depth),
            random.uniform(-max_depth, 0),
            random.uniform(-max_depth, 0),
        ]

        depths = inductiva.generative.procedural.generate_random_map_level(
            x_num, y_num, corner_values, initial_roughness, roughness_factor)

        depths = inductiva.generative.procedural.adjust_map_level(
            depths, percentile_above_water)

        x, y = np.meshgrid(np.linspace(*x_range, x_num),
                           np.linspace(*y_range, y_num),
                           indexing="ij")

        return cls(depths.flatten(), x.flatten(), y.flatten())

    def to_text_file(self, text_file_path: str):
        """Writes the bathymetry to a text file.

        Args:
            text_file_path: Path to the text file.
        """

        np.savetxt(text_file_path, self.depths)

    @property
    def x_range(self) -> Tuple[float]:
        """Returns the range of x values."""

        return (np.min(self.x), np.max(self.x))

    @property
    def y_range(self) -> Tuple[float]:
        """Returns the range of y values."""

        return (np.min(self.y), np.max(self.y))

    def x_uniques(self, sort: bool = False) -> np.ndarray:
        """Returns the unique x values.
        
        Args:
            sort: Whether to sort the unique values.
        """
        x_uniques = np.unique(self.x)
        if sort:
            x_uniques = np.sort(x_uniques)
        return x_uniques

    def y_uniques(self, sort: bool = False) -> np.ndarray:
        """Returns the unique y values.
        
        Args:
            sort: Whether to sort the unique values.
        """
        y_uniques = np.unique(self.y)
        if sort:
            y_uniques = np.sort(y_uniques)
        return y_uniques

    def is_uniform_grid(self) -> bool:
        """Determines whether the bathymetry is defined on a uniform grid."""

        x_uniques_diffs = np.diff(self.x_uniques(sort=True))
        y_uniques_diffs = np.diff(self.y_uniques(sort=True))

        return (np.unique(x_uniques_diffs.round(decimals=2)).size == 1 and
                np.unique(y_uniques_diffs.round(decimals=2)).size == 1)

    def plot(
        self,
        cmap: Optional[str] = None,
        clim: Optional[Tuple[float]] = None,
        path: Optional[str] = None,
        grid_size: int = 500,
        max_distance: float = 20,
    ) -> Union[matplotlib.axes.Axes, None]:
        """Plots the bathymetry.

        The bathymetry is represented as a 2D map of depths, with the x and y
        coordinates of the points where the depths are defined in the axes.

        The bathymetry data is plotted with a color plot on a uniform grid,
        defined by the size `grid_size`.
        
        The data is interpolated from the points where the bathymetry is defined
        to the uniform grid using linear interpolation.

        Points on the uniform grid at a distance to points where the bathymetry
        is defined larger than a threshold distance `max_distance` are not
        plotted.
    
        The plot is produced with matplotlib.

        Args:
            cmap: Colormap to use. Defaults to the matplotlib default colormap.
            clim: Range of depth values to represent in colors. If `None`, the
              range is determined automatically from the minumum and maximum
              depth. Depth values below or above this range are represented with
              the minimum or maximum colors, respectively.
            path: Path to save the plot. If `None`, the plot is not saved, and
              the matplotlib `Axes` object is returned instead.
        """

        # Create uniform grid for interpolation.
        x_grid, y_grid = np.meshgrid(np.linspace(*self.x_range, grid_size),
                                     np.linspace(*self.y_range, grid_size))

        # Interpolate depths to uniform grid.
        interpolator = scipy.interpolate.LinearNDInterpolator(
            (self.x, self.y),
            self.depths,
        )

        depths_grid = interpolator(x_grid, y_grid)

        # Filter out points that are far from bathymetry locations.
        tree = scipy.spatial.KDTree(np.c_[self.x, self.y])
        distance, _ = tree.query(np.c_[x_grid.ravel(), y_grid.ravel()], k=1)
        distance = distance.reshape(x_grid.shape)
        depths_grid[distance > max_distance] = np.nan

        # Plot the bathymetry.
        extent = (
            self.x_range[0],
            self.x_range[1],
            self.y_range[0],
            self.y_range[1],
        )

        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot()

        im = ax.imshow(
            depths_grid,
            cmap=cmap,
            clim=clim,
            origin="lower",
            extent=extent,
        )

        ax.set(
            aspect="equal",
            xlim=self.x_range,
            ylim=self.y_range,
            xlabel="$x$ [m]",
            ylabel="$y$ [m]",
        )

        fig.colorbar(im, ax=ax, label="Depth [m]")

        if path is not None:
            fig.savefig(path)
            matplotlib.pyplot.close(fig)

        else:
            return ax

    def to_uniform_grid(self, x_size: int = 200, y_size: int = 200):
        """Converts the bathymetry to a uniform grid.

        The bathymetry is interpolated to a uniform grid. The grid is defined
        by the ranges of x and y values, and the number of points in each
        direction is configurable via the `x_size` and `y_size` arguments.

        Args:
            x_size: Number of grid points in the x direction.
            y_size: Number of grid points in the y direction.
        """

        x_grid, y_grid = np.meshgrid(np.linspace(*self.x_range, x_size),
                                     np.linspace(*self.y_range, y_size),
                                     indexing="ij")

        depths_grid = scipy.interpolate.griddata(
            (self.x, self.y),
            self.depths,
            (x_grid, y_grid),
            method="linear",
        )

        if np.sum(np.isnan(depths_grid)) > 0:
            raise ValueError(
                "The bathymetry cannot be converted to a uniform grid because "
                "depths are not defined in one or more edge regions of the "
                "domain.")

        return Bathymetry(depths=depths_grid.flatten(),
                          x=x_grid.flatten(),
                          y=y_grid.flatten())
