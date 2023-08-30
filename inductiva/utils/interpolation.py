"""Interpolation utilities."""

from typing import Tuple

import numpy as np
import scipy


def interpolate_to_uniform_grid(x: Tuple[np.ndarray],
                                values: np.ndarray,
                                x_grid: Tuple[np.ndarray],
                                threshold_distance=None):
    """Interpolate values to a uniform grid.

    Args:
        x: Tuple of arrays with the coordinates of the points where values are
          defined.
        values: Array with the values defined at the points defined by `x`.
        x_grid: Tuple of arrays with the coordinates of the points where values
          will be interpolated.
        threshold_distance: If not None, points in `x_grid` that are further
          than `threshold_distance` from the closest point in `x` will be set to
          NaN.
    
    Returns:
        Array with the interpolated values.
    """

    interpolator = scipy.interpolate.LinearNDInterpolator(x, values)

    values_grid = interpolator(x_grid)

    if threshold_distance is not None:
        # Filter out points that are far from x.
        tree = scipy.spatial.KDTree(np.transpose(x))

        # Obtain distance between each point on the uniform grid and the closest
        # location where values are defined.
        x_grid_shape = np.shape(x_grid)

        distance, _ = tree.query(
            np.transpose(np.reshape(x_grid, (x_grid_shape[0], -1))),
            k=1,
        )
        distance = distance.reshape(x_grid_shape[1:])

        # Set values to NaN for points that are far from locations where values
        # are defined.
        values_grid[distance > threshold_distance] = np.nan

    return values_grid
