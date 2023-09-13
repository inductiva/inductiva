"""Interpolation utilities."""

from typing import Literal, Optional, Tuple

import numpy as np
try:
    import scipy
except ImportError:
    scipy = None


def interpolate_to_coordinates(
    coordinates: Tuple[np.ndarray],
    values: np.ndarray,
    interpolation_coordinates: Tuple[np.ndarray],
    method: Literal["nearest", "linear"] = "linear",
    threshold_distance: Optional[float] = None,
):
    """Interpolate values defined at coordinates x to interpolation coordinates.

    Args:
        coordinates: Tuple of arrays with the coordinates of the points where
          values are defined.
        values: Array with the values defined at the points defined by `x`.
        interpolation_coordinates: Tuple of arrays with the coordinates of the
          points where values are interpolated to.
        threshold_distance: If not `None`, points in `interpolation_coordinates`
          that are farther than `threshold_distance` from the closest point in
          `coordinates` are set to `np.NAN`.
    
    Returns:
        Array with the interpolated values. Its shape matches that of the
        interpolation coordinate arrays.
    """

    if method == "nearest":
        interpolator = scipy.interpolate.NearestNDInterpolator(
            coordinates, values)
    elif method == "linear":
        interpolator = scipy.interpolate.LinearNDInterpolator(
            coordinates, values)
    else:
        raise ValueError("Invalid interpolation method.")

    values_interpolation = interpolator(interpolation_coordinates)

    if threshold_distance is not None:
        # Filter out points that are far from x.
        tree = scipy.spatial.KDTree(np.transpose(coordinates))

        # Obtain distance between each interpolation point and the closest
        # location where values are defined.
        interpolation_coordinates_shape = np.shape(interpolation_coordinates)

        # Queries to the KDTree require that interpolation coordinates are
        # provided in shape num_points x dim_coordinates, where num_points
        # is the number of points to query and dim_coordinates is the number
        # of coordinates per point.
        # Reshape interpolation coordinates to match that shape.
        interpolation_coordinates_reshaped = np.reshape(
            interpolation_coordinates, (interpolation_coordinates_shape[0], -1))
        interpolation_coordinates_reshaped = np.transpose(
            interpolation_coordinates_reshaped)
        distance, _ = tree.query(interpolation_coordinates_reshaped, k=1)
        distance = distance.reshape(interpolation_coordinates_shape[1:])

        # Set values to NaN for points that are far from locations where values
        # are defined.
        values_interpolation[distance > threshold_distance] = np.nan

    return values_interpolation
