"""Interpolation utilities."""

from typing import Literal, Optional, Tuple

import numpy as np
import scipy


def interpolate_to_coordinates(
    x: Tuple[np.ndarray],
    values: np.ndarray,
    x_interpolation: Tuple[np.ndarray],
    method: Literal["nearest", "linear"] = "linear",
    threshold_distance: Optional[float] = None,
):
    """Interpolate values defined at coordinates x to interpolation coordinates.

    Args:
        x: Tuple of arrays with the coordinates of the points where values are
          defined.
        values: Array with the values defined at the points defined by `x`.
        x_interpolation: Tuple of arrays with the coordinates of the points
          where values will be interpolated.
        threshold_distance: If not None, points in `x_interpolation` that are
          further than `threshold_distance` from the closest point in `x` are
          set to NaN.
    
    Returns:
        Array with the interpolated values. Its shape matches that of the
        interpolation coordinate arrays.
    """

    if method == "nearest":
        interpolator = scipy.interpolate.NearestNDInterpolator(x, values)
    elif method == "linear":
        interpolator = scipy.interpolate.LinearNDInterpolator(x, values)
    else:
        raise ValueError("Invalid interpolation method.")

    values_interpolation = interpolator(x_interpolation)

    if threshold_distance is not None:
        # Filter out points that are far from x.
        tree = scipy.spatial.KDTree(np.transpose(x))

        # Obtain distance between each interpolation point and the closest
        # location where values are defined.
        x_interpolation_shape = np.shape(x_interpolation)

        distance, _ = tree.query(
            np.transpose(
                np.reshape(x_interpolation, (x_interpolation_shape[0], -1))),
            k=1,
        )
        distance = distance.reshape(x_interpolation_shape[1:])

        # Set values to NaN for points that are far from locations where values
        # are defined.
        values_interpolation[distance > threshold_distance] = np.nan

    return values_interpolation
