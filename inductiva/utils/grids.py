"""Utility functions to manipulate grid data."""
import typing

import numpy as np
try:
    import scipy
except ImportError:
    scipy = None


def get_meshgrid(x_range: typing.Sequence[float],
                 y_range: typing.Sequence[float], x_num: int, y_num: int):

    x_grid = np.linspace(*x_range, num=x_num)
    y_grid = np.linspace(*y_range, num=y_num)

    x_grid, y_grid = np.meshgrid(x_grid, y_grid, indexing="ij")

    return x_grid, y_grid


def reshape_map(x_num: int, y_num: int, map_level: np.ndarray):
    """Reshape a map into a different resolution.
    
    Args:
        x_num: New number of grid points in x-direction.
        y_num: New number of grid points in y-direction.
        map_level: Values on the z-axis of the previous grid.
    """

    x_grid_prev, y_grid_prev = np.meshgrid(np.arange(map_level.shape[0]),
                                           np.arange(map_level.shape[1]),
                                           indexing="ij")

    x_grid_new, y_grid_new = np.meshgrid(np.arange(x_num),
                                         np.arange(y_num),
                                         indexing="ij")

    return scipy.interpolate.griddata(
        (x_grid_prev.flatten(), y_grid_prev.flatten()),
        map_level.flatten(),
        (x_grid_new, y_grid_new),
        method="linear",
    )
