"""Bathymetry class and utilities."""

import random
from typing import Optional, Sequence, Tuple, Union

from absl import logging

import numpy as np
try:
    import matplotlib
    import utm
except ImportError:
    matplotlib = None
    utm = None

import inductiva
from inductiva.utils import optional_deps


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

    @optional_deps.needs_coastal_extra_deps
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
    @optional_deps.needs_coastal_extra_deps
    def from_bot_file(
        cls,
        bot_file_path: str,
        x_range: Sequence[float],
        y_range: Sequence[float],
    ):
        """Creates a `Bathymetry` object from a bot file.

        The depth values are read from a bot file, i.e. a text file with a 2D
        array with the depths, in meters. The first and second dimensions of the
        array in the text file (i.e. rows and columns) correspond to the x and y
        directions, respectively.

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
    @optional_deps.needs_coastal_extra_deps
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
    @optional_deps.needs_coastal_extra_deps
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
        random_seed=None,
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
            random_seed: Random seed to use.
        """

        random.seed(random_seed)

        corner_values = [
            random.uniform(0, max_depth),
            random.uniform(0, max_depth),
            random.uniform(-max_depth, 0),
            random.uniform(-max_depth, 0),
        ]

        depths = inductiva.generative.procedural.generate_random_map_level(
            x_num,
            y_num,
            corner_values,
            initial_roughness,
            roughness_factor,
            random_seed=random_seed,
        )

        depths = inductiva.generative.procedural.adjust_map_level(
            depths, percentile_above_water)

        x, y = np.meshgrid(np.linspace(*x_range, x_num),
                           np.linspace(*y_range, y_num),
                           indexing="ij")

        return cls(depths.flatten(), x.flatten(), y.flatten())

    @optional_deps.needs_coastal_extra_deps
    def to_bot_file(self, bot_file_path: str, depths_grid: np.ndarray):
        """Writes the bathymetry to a bot file.

        The depth values are interpolated to a regular grid and written to a bot
        file. The grid size is determined automatically from the range and the
        number of unique x and y values.

        Args:
            text_file_path: Path to the text file.
        """

        np.savetxt(bot_file_path, depths_grid)

    @property
    @optional_deps.needs_coastal_extra_deps
    def x_range(self) -> Tuple[float]:
        """Returns the range of x values."""

        return (np.min(self.x), np.max(self.x))

    @property
    @optional_deps.needs_coastal_extra_deps
    def y_range(self) -> Tuple[float]:
        """Returns the range of y values."""

        return (np.min(self.y), np.max(self.y))

    @optional_deps.needs_coastal_extra_deps
    def x_ptp(self) -> float:
        """Returns the peak-to-peak range (max - min) of x values."""

        return np.ptp(self.x)

    @optional_deps.needs_coastal_extra_deps
    def y_ptp(self) -> float:
        """Returns the peak-to-peak range (max - min) of y values."""

        return np.ptp(self.y)

    @optional_deps.needs_coastal_extra_deps
    def x_num(self) -> int:
        """Returns the length of unique x values."""

        return int(len(self.x_uniques()))

    @optional_deps.needs_coastal_extra_deps
    def y_num(self) -> int:
        """Returns the length of unique y values."""

        return int(len(self.y_uniques()))

    @optional_deps.needs_coastal_extra_deps
    def x_uniques(self, sort: bool = False) -> np.ndarray:
        """Returns the unique x values.

        Args:
            sort: Whether to sort the unique values.
        """
        x_uniques = np.unique(self.x)
        if sort:
            x_uniques = np.sort(x_uniques)
        return x_uniques

    @optional_deps.needs_coastal_extra_deps
    def y_uniques(self, sort: bool = False) -> np.ndarray:
        """Returns the unique y values.

        Args:
            sort: Whether to sort the unique values.
        """
        y_uniques = np.unique(self.y)
        if sort:
            y_uniques = np.sort(y_uniques)
        return y_uniques

    @optional_deps.needs_coastal_extra_deps
    def depths_grid(self) -> np.ndarray:
        """Returns the depths as a 2D grid."""

        return self.depths.reshape((self.x_num(), self.y_num()))

    @optional_deps.needs_coastal_extra_deps
    def crop(self,
             x_range: Sequence[float],
             y_range: Sequence[float],
             remove_offset=True):
        """Crops the bathymetry to a given range of x and y values.

        Args:
            x_range: The range of x values, in meters.
            y_range: The range of y values, in meters.
        """
        x_min, x_max = x_range
        y_min, y_max = y_range

        mask = (self.x >= x_min) & (self.x <= x_max) & \
               (self.y >= y_min) & (self.y <= y_max)

        if np.sum(mask) == 0:
            raise ValueError(
                "The bathymetry cannot be cropped because no points are "
                "defined in the given coordinate ranges.")

        x = self.x[mask]
        y = self.y[mask]
        depths = self.depths[mask]

        if remove_offset:
            x = x - x.min()
            y = y - y.min()

        return Bathymetry(
            depths=depths,
            x=x,
            y=y,
        )

    @optional_deps.needs_coastal_extra_deps
    def is_uniform_grid(self) -> bool:
        """Determines whether the bathymetry is defined on a uniform grid."""

        x_uniques_diffs = np.diff(self.x_uniques(sort=True))
        y_uniques_diffs = np.diff(self.y_uniques(sort=True))

        return (np.unique(x_uniques_diffs.round(decimals=2)).size == 1 and
                np.unique(y_uniques_diffs.round(decimals=2)).size == 1)

    @optional_deps.needs_coastal_extra_deps
    def plot(
        self,
        show: bool = False,
        cmap: str = "coolwarm",
        clim: Optional[Tuple[float]] = None,
        path: Optional[str] = None,
        x_resolution: float = 10,
        y_resolution: float = 10,
        resize: bool = False,
        threshold_distance: float = 20,
    ) -> Union["matplotlib.axes.Axes", None]:
        """Plots the bathymetry.

        The bathymetry is represented as a 2D map of depths, with the x and y
        coordinates of the points where the depths are defined in the axes. The
        plot shows the height of the various points in the bathymetry is
        defined, with the opposite signal convention of the bathymetry (i.e,
        the bellow water level points are given by negative height instead of
        positive depth).

        The bathymetry data is plotted with a color plot on a grid with uniform
        spacing in the x and y directions. If the user does not specify the
        resolution, the plotted mesh has the bathymetry resolution. If
        necessary, the user can control the spacing in the x and y directions
        using the `x_resolution` and `y_resolution` arguments.

        For non-uniform grids or if the resoltuion is specified, the data is
        interpolated from the points where the bathymetry is defined
        to the resized uniform grid using linear interpolation.

        Points on the uniform grid at a distance larger than a threshold
        distance `threshold_distance` from points where the bathymetry is
        defined are omitted.

        The plot is produced with matplotlib.

        Args:
            cmap: Colormap to use. Defaults to the matplotlib default colormap.
            clim: Range of depth values to represent in colors. If `None`, the
              range is determined automatically from the minumum and maximum
              depth. Depth values below or above this range are represented with
              the minimum or maximum colors, respectively.
            path: Path to save the plot. If `None`, the plot is not saved, and
              the matplotlib `Axes` object is returned instead.
            x_resolution: Resolution, in meters, of the plotting grid in the x
              direction. Only used if the bathymetry is not uniform or if we
              are resizing it.
            y_resolution: Resolution, in meters, of the plotting grid in the y
              direction. Only used if the bathymetry is not uniform or if we
              are resizing it.
            resize: Whether to resize the bathymetry before plotting it. When
                `False`, the bathymetry is plotted on the original grid if
                it is already uniform.
            threshold_distance: Threshold distance to filter out points on the
              uniform grid that are far from points where the bathymetry is
              defined.
        """

        if self.is_uniform_grid() and not resize:
            x_size = self.x_num()
            y_size = self.y_num()

            depths_grid = self.depths.reshape((x_size, y_size))

        else:
            logging.info("Plotting the bathymetry on resized uniform grid...")

            depths_grid, _ = self._interpolate_to_uniform_grid(
                x_resolution,
                y_resolution,
                threshold_distance=threshold_distance,
                nullable=True,
            )

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
            -1 * depths_grid.transpose(),
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

        fig.colorbar(im, ax=ax, label="Height [m]")

        if path is not None:
            fig.savefig(path)
            matplotlib.pyplot.close(fig)

        else:
            if show:
                matplotlib.pyplot.show()
            return ax

    @optional_deps.needs_coastal_extra_deps
    def to_uniform_grid(
        self,
        x_resolution: float = 2,
        y_resolution: float = 2,
        fill_value: Optional[Union[float, str]] = None,
    ):
        """Converts the bathymetry to a uniform grid.

        The bathymetry is interpolated to a grid with uniform resolution (i.e.
        spacing) in the x and y directions. Linear interpolation is used to
        interpolate the bathymetry to the grid.

        Grid points for which no interpolation is possible may be filled with
        a constant or with the nearest depth value.

        Args:
            x_resolution: Resolution, in meters, of the grid in the x direction.
            y_resolution: Resolution, in meters, of the grid in the y direction.
            fill_value: Value to fill the grid points for which no interpolation
              is possible. If "nearest", the nearest depth value is used.
        """

        depths_grid, (x_grid, y_grid) = self._interpolate_to_uniform_grid(
            x_resolution, y_resolution, fill_value=fill_value, nullable=False)

        return Bathymetry(depths=depths_grid.flatten(),
                          x=x_grid.flatten(),
                          y=y_grid.flatten())

    @optional_deps.needs_coastal_extra_deps
    def _interpolate_to_uniform_grid(
        self,
        x_resolution: float,
        y_resolution: float,
        threshold_distance: Optional[float] = None,
        fill_value: Optional[Union[float, str]] = None,
        nullable: bool = False,
    ):
        """Interpolates the bathymetry to a uniform grid.

        Args:
            x_resolution: Resolution, in meters, of the grid in the x direction.
            y_resolution: Resolution, in meters, of the grid in the y direction.
            threshold_distance: Threshold distance to filter out points on the
              uniform grid that are far from points where the bathymetry is
              defined.
            fill_value: Value to fill the grid points for which no interpolation
              is possible. If "nearest", the nearest depth value is used.
            nullable: Whether to allow the bathymetry to be undefined in some
              grid points. If `False`, an error is raised if the bathymetry is
              undefined in one or more grid points.

        Returns:
            depths_grid: A 2D array with the depths on the uniform grid.
            (x_grid, y_grid): The x and y coordinates of the points where the
              depths are defined, respectively.
        """
        # Determine grid size based on ranges and resolution.
        x_size = int(self.x_ptp() / x_resolution)
        y_size = int(self.y_ptp() / y_resolution)

        logging.info(
            "Interpolating the bathymetry to a grid of dimensions:\n"
            "- grid resolution %.2f x %.2f (m x m) \n"
            "- grid size %d x %d", x_resolution, y_resolution, x_size, y_size)

        # Create uniform grid for interpolation.
        (x_grid, y_grid) = inductiva.utils.grids.get_meshgrid(
            x_range=self.x_range,
            y_range=self.y_range,
            x_num=x_size,
            y_num=y_size,
        )

        depths_grid = inductiva.utils.interpolation.interpolate_to_coordinates(
            coordinates=(self.x, self.y),
            values=self.depths,
            interpolation_coordinates=(x_grid, y_grid),
            method="linear",
            threshold_distance=threshold_distance)

        nan_mask = np.isnan(depths_grid)

        if fill_value is not None:
            if fill_value == "nearest":
                nearest_depths_grid = \
                    inductiva.utils.interpolation.interpolate_to_coordinates(
                        coordinates=(self.x, self.y),
                        values=self.depths,
                        interpolation_coordinates=(x_grid, y_grid),
                        method="nearest")

                depths_grid[nan_mask] = nearest_depths_grid[nan_mask]
            else:
                depths_grid[nan_mask] = fill_value

        nan_mask = np.isnan(depths_grid)

        if not nullable and np.sum(nan_mask) > 0:
            raise ValueError(
                "The bathymetry cannot be interpolated to a uniform grid "
                "because depths are not defined in one or more grid points.")

        return depths_grid, (x_grid, y_grid)
