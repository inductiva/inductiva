"""Bathymetry class and utilities."""

import random
from typing import Optional, Sequence, Tuple, Union

import matplotlib
import numpy as np

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
    ) -> Union[matplotlib.axes.Axes, None]:
        """Plots the bathymetry.

        The bathymetry is represented as a 2D map of depths, with the x and y
        coordinates of the points where the depths are defined in the axes.

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

        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot()

        pc = ax.tripcolor(
            self.x,
            self.y,
            self.depths,
            cmap=cmap,
            clim=clim,
        )
        ax.set(
            aspect="equal",
            xlim=self.x_range,
            ylim=self.y_range,
            xlabel="$x$ [m]",
            ylabel="$y$ [m]",
        )

        fig.colorbar(pc, ax=ax, label="Depth [m]")

        if path is not None:
            fig.savefig(path)
            matplotlib.pyplot.close(fig)

        else:
            return ax
