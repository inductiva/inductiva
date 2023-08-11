"""Utility functions to manipulate grid data."""
import typing

import scipy
import numpy as np


def get_meshgrid(x_range: typing.Sequence[float],
                 y_range: typing.Sequence[float], num_list: typing.List[int]):
    """Create a grid of x and y values.

    Given a [x_range, y_range] and a resolution for
    each range, we return a meshgrid of x and y values.
    That is, x_grid[i, j] = x_range[i] and
    y_grid[i, j] = y_range[j].

    Args:
        x_range: The range of x values, in meters.
        y_range: The range of y values, in meters.
        num_list: Number of grid points in both
            directions [x_num, y_num]
    """

    x_grid = np.linspace(*x_range, num=num_list[0])
    y_grid = np.linspace(*y_range, num=num_list[1])

    x_grid, y_grid = np.meshgrid(x_grid, y_grid, indexing="ij")

    return x_grid, y_grid


def interpolate_between_grids(x_range: typing.Sequence[float],
                              y_range: typing.Sequence[float],
                              prev_num_list: typing.List[int],
                              new_num_list: typing.List[int],
                              z_array: np.ndarray):
    """Interpolate between two grid with different resolution.
    
    Args:
        x_range: The range of x values.
        y_range: The range of y values.
        prev_num_list: Number of grid points in both directions of the
            previous grid [x_num, y_num].
        new_num_list: Number of grid points in both directions of the
            new grid [x_num, y_num].
        z_array: Values on the z-axis of the previous grid.
    """

    x_grid_prev, y_grid_prev = get_meshgrid(x_range=x_range,
                                            y_range=y_range,
                                            num_list=prev_num_list)

    x_grid_new, y_grid_new = get_meshgrid(x_range=x_range,
                                          y_range=y_range,
                                          num_list=new_num_list)

    new_z_array = scipy.interpolate.griddata(
        (x_grid_prev.flatten(), y_grid_prev.flatten()),
        z_array.flatten(),
        (x_grid_new, y_grid_new),
        method="linear",
    )

    return new_z_array
