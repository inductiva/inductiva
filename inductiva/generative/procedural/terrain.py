"""Functions to generate elevation data for a general grid."""

import math
import typing

import numpy as np

import inductiva
from inductiva.generative import procedural


def generate_random_terrain(
    x_range: typing.Sequence[float],
    y_range: typing.Sequence[float],
    num_list: typing.List[int],
    corner_values: typing.Sequence[float],
    initial_roughness: float = 1,
    roughness_factor: float = 0.5,
    percentile_translate_terrain: float = 0,
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
        x_range: The range of x values, in meters.
        y_range: The range of y values, in meters.
        num_list: The number of points in the x and y directions.
        corner_values: Sequence of 4 values establishing the elevation of the 
            grid corners. The order refers to top-left, top-right,
            bottom-left, and bottom-right.
        initial_roughness: Initial roughness value, in meters. Controls the
            initial range of randomness of the Diamond-Square algorithm.
        roughness_factor: Roughness factor. Must be between 0 and 1.
            Controls the rate at which the range of randomness of the
            Diamond-Square algorithm decreases over iterations.
        percentile_translate_terrain: Percentile of the elevation that must
            be above xy-plane. Must be between -100 and 100, where negative
            values set the terrain lower and positive above.
    """

    # Determine the minimum n such that 2^n + 1 >= max(x_num, y_num).
    n_power = int(math.log2(max(num_list) - 1)) + 1

    size_square = 2**n_power + 1

    # Create elevation for a square grid with side resolution=size_square.
    z_elevation = procedural.diamond_square.create_random_array(
        size=size_square,
        corner_values=corner_values,
        initial_roughness=initial_roughness,
        roughness_factor=roughness_factor,
        random_seed=random_seed)

    # Adjust terrain to ensure that a given percentage of the terrain elevation
    # is above or below the xy-plane.
    percentile_translate = np.percentile(z_elevation,
                                         abs(percentile_translate_terrain))
    z_elevation += np.sign(percentile_translate_terrain) * percentile_translate

    z_elevation = inductiva.utils.grids.interpolate_between_grids(
        x_range=x_range,
        y_range=y_range,
        prev_num_list=[size_square, size_square],
        new_num_list=num_list,
        z_array=z_elevation)

    return z_elevation
