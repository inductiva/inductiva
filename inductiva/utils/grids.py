"""Utility functions to manipulate grid data."""
import typing

import scipy
import numpy as np


def get_meshgrid(x_range: typing.Sequence[float],
                 y_range: typing.Sequence[float], x_num: int, y_num: int):

    x_grid = np.linspace(*x_range, num=x_num)
    y_grid = np.linspace(*y_range, num=y_num)

    x_grid, y_grid = np.meshgrid(x_grid, y_grid, indexing="ij")

    return x_grid, y_grid


def interpolate_between_grids(x_num: int, y_num: int, z_array: np.ndarray):
    """Interpolate between two grid with different resolution.
    
    Args:
        x_num: New number of grid points in x-direction.
        y_num: New number of grid points in y-direction.
        z_array: Values on the z-axis of the previous grid.
    """

    x_grid_prev, y_grid_prev = np.meshgrid(np.arange(z_array.shape[0]),
                                           np.arange(z_array.shape[1]),
                                           indexing="ij")

    x_grid_new, y_grid_new = np.meshgrid(np.arange(x_num),
                                         np.arange(y_num),
                                         indexing="ij")

    new_z_array = scipy.interpolate.griddata(
        (x_grid_prev.flatten(), y_grid_prev.flatten()),
        z_array.flatten(),
        (x_grid_new, y_grid_new),
        method="linear",
    )

    return new_z_array
