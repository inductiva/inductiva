"""Functions to generate elevation data for a general grid."""

import math
import typing

import numpy as np
import scipy

import inductiva


def generate_grid_elevation(
    x_range: typing.Sequence[float],
    y_range: typing.Sequence[float],
    x_num: int,
    y_num: int,
    corner_values: typing.Sequence[float],
    initial_roughness: float = 1,
    roughness_factor: float = 0.5,
    percentile_translate_terrain: float = 0,
):
    """Generate a set of elevation values for xy-grid.
    
    This function generalises the diamond-square algorithm for grids that are
    not square and multiple of 2^n+1.

    A grid of random elevations of shape `(x_num, y_num)` is generated as
    follows:
    1. A square grid of random elevations is generated using the Diamond-Square
    algorithm. The randomness of the algorithm is controlled
    by the `initial_roughness` and `roughness_factor` parameters, that set
    the initial range of randomness and the rate at which it decreases over
    iterations of the algorithm, respectively. The size of the resultnig
    grid is `n`, where `n` is the smallest integer that satisfies
    `2^n + 1 >= max(x_num, y_num)`.
    2. The elevations of the square grid are interpolated to a grid of shape
    `(x_num, y_num)`.

    Args:
        x_range: The range of x values, in meters.
        y_range: The range of y values, in meters.
        x_num: Number of grid points in the x direction.
        y_num: Number of grid points in the y direction.
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
    n_power = int(math.log2(max(x_num, y_num) - 1)) + 1

    size_square = 2**n_power + 1

    # Create elevation for a square grid with side resolution=size_square.
    z_elevation = inductiva.generative.diamond_square.create_random_array(
        size=size_square,
        corner_values=corner_values,
        initial_roughness=initial_roughness,
        roughness_factor=roughness_factor)

    # Adjust terrain to ensure that a given percentage of the terrain elevation
    # is above or below the xy-plane.
    percentile_translate = np.percentile(z_elevation,
                                         abs(percentile_translate_terrain))
    z_elevation += np.sign(percentile_translate_terrain) * percentile_translate

    # Interpolate from the square grid to the desired grid.
    x_square = np.linspace(*x_range, size_square)
    y_square = np.linspace(*y_range, size_square)

    x_square, y_square = np.meshgrid(x_square, y_square, indexing="ij")

    x_grid = np.linspace(*x_range, x_num)
    y_grid = np.linspace(*y_range, y_num)

    x_grid, y_grid = np.meshgrid(x_grid, y_grid, indexing="ij")

    z_elevation = scipy.interpolate.griddata(
        (x_square.flatten(), y_square.flatten()),
        z_elevation.flatten(),
        (x_grid, y_grid),
        method="linear",
    )

    return x_grid, y_grid, z_elevation
