"""Functions to generate elevation data for a general grid."""

import math
import typing

import numpy as np

import inductiva
from inductiva.generative import procedural


def generate_random_map_level(
    x_num: int,
    y_num: int,
    corner_values: typing.Sequence[float],
    initial_roughness: float = 1,
    roughness_factor: float = 0.5,
    random_seed: int = None,
):
    """Generate a set of elevation values for xy-grid.
    
    This function generalises the diamond-square algorithm for grids that are
    not square and multiple of 2^n+1.

    A grid of random elevations of shape `(x_num, y_num)` is generated as
    follows:
    1. A square grid of random elevations is generated using the Diamond-Square
    algorithm. The size of the resultnig grid is `n`, where `n` is
    the smallest integer that satisfies `2^n + 1 >= max(x_num, y_num)`.
    2. The elevations of the square grid are interpolated to a grid of shape
    `(x_num, y_num)`.

    Args:
        x_num: The number of points in the x direction.
        y_num: The number of points in the y direction.
        corner_values: Sequence of 4 values establishing the elevation of the 
            grid corners. The order refers to top-left, top-right,
            bottom-left, and bottom-right.
        initial_roughness: Initial roughness value, in meters. Controls the
            initial range of randomness of the Diamond-Square algorithm.
        roughness_factor: Roughness factor. Must be between 0 and 1.
            Controls the rate at which the range of randomness of the
            Diamond-Square algorithm decreases over iterations.
        random_seed: Random seed to use for the Diamond-Square algorithm.
    """

    # Determine the minimum n such that 2^n + 1 >= max(x_num, y_num).
    n_power = int(math.log2(max(x_num, y_num) - 1)) + 1

    size_square = 2**n_power + 1

    # Create elevation for a square grid with side resolution=size_square.
    map_level = procedural.diamond_square.create_random_array(
        size=size_square,
        corner_values=corner_values,
        initial_roughness=initial_roughness,
        roughness_factor=roughness_factor,
        random_seed=random_seed)

    return inductiva.utils.grids.reshape_map(x_num=x_num,
                                             y_num=y_num,
                                             map_level=map_level)


def adjust_map_level(map_level: np.ndarray, percentile_translate_map: float):
    """Adjust a map level based on percentile.
    
    Args:
        map_level: The map level to adjust.
        percentile_translate_map: The percentile of the map level to adjust.
        """
    # Adjust a map_level to ensure that a given percentage of it
    # is above or below the xy-plane.
    percentile_translate = np.percentile(map_level,
                                         abs(percentile_translate_map))
    map_level -= np.sign(percentile_translate_map) * percentile_translate

    return map_level
