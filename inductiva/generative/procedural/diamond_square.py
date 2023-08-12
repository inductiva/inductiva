"""Utilities to generate random 2D arrays with the Diamond-Square algorithm."""

import math
import random
from typing import Sequence

from absl import logging

import numpy as np


def kernel_avg(kernel_offsets,
               random_array,
               center_x,
               center_y,
               kernel_scale,
               enforce_inner_indexes=True):
    """Averaging kernel, given list of relative offsets and scaling factor."""

    # Scale to the distance and then apply offset
    scaled_kernel_offsets = np.multiply(kernel_offsets, kernel_scale)
    src_idxs = scaled_kernel_offsets + [center_x, center_y]

    # Make sure we are only getting indexes inside the map.
    if enforce_inner_indexes:
        all_positive = (src_idxs[:, 0] >= 0) & (src_idxs[:, 1] >= 0)
        max_x, max_y = random_array.shape
        all_bounded = (src_idxs[:, 0] < max_x) & (src_idxs[:, 1] < max_y)
        # Filter out to make sure all indices are inside the map.
        src_idxs = src_idxs[all_positive & all_bounded]

    # Fetch elements by list of x and list of y coords and return the mean.
    elements = random_array[src_idxs[:, 0], src_idxs[:, 1]]
    return np.mean(elements)


def diamond_step(random_array, inter_center_distance, noise_range, random_seed):
    """The diamond step performed at a given inter center distance on a map.
                            +-----------+
                            | D |   | D |
                            |---+---+---|
                            |   | X |   |
                            |---+---+---|
                            | D |   | D |
                            +-----------+
    For each square in the array, set the midpoint to be the average of the four
    corner points plus a random value."""

    diamond_offsets = [[-1, -1], [-1, 1], [1, 1], [1, -1]]
    map_length = random_array.shape[0]
    kernel_scale = inter_center_distance // 2
    random.seed(random_seed)

    for i in range(kernel_scale, map_length, inter_center_distance):
        for j in range(kernel_scale, map_length, inter_center_distance):
            random_array[i, j] = kernel_avg(diamond_offsets, random_array, i, j,
                                            kernel_scale)
            random_array[i, j] += random.uniform(-noise_range, noise_range)


def square_step(random_array, inter_center_distance, noise_range, random_seed):
    """The square step performed at a given inter center distance on a map.
                            +-----------+
                            |   | S |   |
                            |---+---+---|
                            | S | X | S |
                            |---+---+---|
                            |   | S |   |
                            +-----------+
    For each diamond in the array, set the midpoint to be the average of the
    four corner points plus a random value."""

    square_offsets = [[-1, 0], [0, -1], [1, 0], [0, 1]]
    map_length = random_array.shape[0]
    kernel_scale = inter_center_distance // 2
    random.seed(random_seed)
    # Square Step, rows
    for i in range(kernel_scale, map_length, inter_center_distance):
        for j in range(0, map_length, inter_center_distance):
            random_array[i, j] = kernel_avg(square_offsets, random_array, i, j,
                                            kernel_scale)
            random_array[i, j] += random.uniform(-noise_range, noise_range)

    # Square Step, cols
    for i in range(0, map_length, inter_center_distance):
        for j in range(kernel_scale, map_length, inter_center_distance):
            random_array[i, j] = kernel_avg(square_offsets, random_array, i, j,
                                            kernel_scale)
            random_array[i, j] += random.uniform(-noise_range, noise_range)


def iterate_diamond_square(initial_condition: np.ndarray,
                           initial_roughness: float = 1,
                           roughness_factor: float = 0.5,
                           random_seed: int = None):
    """Performs the diamond-square iteration."""

    if initial_condition.shape[0] != initial_condition.shape[1]:
        logging.error(
            'The initial array must be a square. The current shape is: %f x %f',
            initial_condition.shape[0], initial_condition.shape[1])

    logging.info('Generating random array with the Diamond-Square Algorithm...')

    map_length = initial_condition.shape[0]
    inter_center_distance = map_length - 1
    roughness = initial_roughness
    random_array = initial_condition

    while inter_center_distance > 1:

        diamond_step(random_array, inter_center_distance, roughness,
                     random_seed)
        square_step(random_array, inter_center_distance, roughness, random_seed)

        # Go to a finer grid, and decrease the amount of noise
        inter_center_distance //= 2
        roughness *= roughness_factor

    return random_array


def create_initial_condition(size: int,
                             corner_values: Sequence[float]) -> np.ndarray:
    """Creates a square array with the given corner values."""

    initial_condition = np.zeros((size, size), dtype=np.float64)
    corner_idxs = [(0, 0), (0, -1), (-1, 0), (-1, -1)]

    for idx, value in zip(corner_idxs, corner_values):
        initial_condition[idx] = value

    return initial_condition


def create_random_array(size: int,
                        corner_values: Sequence[float],
                        initial_roughness: float = 1,
                        roughness_factor: float = 0.5,
                        random_seed: int = None) -> np.ndarray:
    """Creates a random square array.

    The random array is generated with an implementation of the Diamond-Square
    Algorithm.

    The description of the Diamond-Square algorithms can be found at
    https://en.wikipedia.org/wiki/Diamond-square_algorithm.

    We start with a n x n map, (with n = 2^k + 1) whose 4 corners have certain
    values. Thus, we first apply a diamond step, to fill in the center of that
    map:
        A - - - - - - - - - B                         A - - - - - - - - - B  ---
        - - - - - - - - - - -                         -\\ - - - - - - - //-   ^
        - - - - - - - - - - -                         - -\\ - - - - - //- -   |
        - - - - - - - - - - -                         - - -\\ - - - //- - -   |
        - - - - - - - - - - -                         - - - -\\ - //- - - -   v
        - - - - - - - - - - -   -- diamond step -->   - - - - - X - - - - -  ---
        - - - - - - - - - - -                         - - - -// -\\ - - - -
        - - - - - - - - - - -                         - - -// - - -\\ - - -
        - - - - - - - - - - -                         - -// - - - - -\\ - -
        - - - - - - - - - - -                         -// - - - - - - -\\ -
        C - - - - - - - - - D                         C - - - - - - - - - D
                                                      |<------->|
                                                      kernel_scale
    Then:
        * X = AVG(A, B, C, D)
    AVG is an "averaging kernel". This kernel has a property that we will name
    as kernel_scale. The kernel_scale is the distance (in rows/columns) between
    the point to be computed and the sources (from a diamond or a square).

    Then, we perform a square step. Observe (below) that there are 4 points,
    N, S, W and E we need to fill at the edges of the map (so we can have more
    squares):
                            A - - - > N < - - - B
                            | - - - - ^ - - - - |
                            | - - - - | - - - - |
                            | - - - - | - - - - |
                            v - - - - | - - - - v
                            W < - - - X - - - > E
                            ^ - - - - | - - - - ^
                            | - - - - | - - - - |
                            | - - - - | - - - - |
                            | - - - - v - - - - |
                            C - - - > S < - - - D
    These are computed by considering first the rows (for N and S) and then the
    columns (for W and E). These will be the computations:
        * N = AVG(A, B, X)
        * S = AVG(C, D, X)
        * W = AVG(A, C, X)
        * E = AVG(A, B, X)
    It's important to introduce another concept: inter center distance. This is
    the distance, in number of rows or number of columns, between the points
    that get averaged. For example, in the iteration we have just showed the
    inter center distance the width/height of the map (because we are at the
    first iteration).

    In general:
        inter_center_distance = 2 * kernel_scale + 1

    Conversely:
        kernel_scale = inter_center_distance // 2

    The diamond-square iteration is applied multiple times with increasing
    degrees of detail, which will end up covering the entire map. For that, the
    inter center distance is halved at every new diamond-square iteration until
    we get to a distance of 1. At that point, all elements of the map grid
    have been generated.

    This averaging procedure has a natural smoothing effect that ensures
    continuity in the terrain. There is, however, some randomness being added at
    every moment, to ensure we don't get a simple smooth profile. However, the
    amount of noise introduced at every step of the iteration allow us to get
    both realistic randomness and, at the same time, realistic continuity.

    Args:
        size: Array size. Must be 2^k + 1. The shape of the array is (size,
          size).
        initial_roughness: Roughness applied in the first iteration of the
          Diamond-Square algorithm. The roughness is a measure of the amount of
          randomness that is added to the array in each step of the algorithm.
        roughness_factor: Fractional value in range (0, 1) by which the initial
          roughness is multiplied at every iteration of the Diamond-Square
          algorithm. Determines how smooth the final array is. Smaller values
          lead to smoother arrays.
        random_seed: Random seed used to generate the random array.
    Returns:
        Random array.
    """

    if not math.log2(size - 1).is_integer():
        raise ValueError('`size` must be 2^k + 1.')

    if len(corner_values) != 4:
        raise ValueError('`corner_values` must be a sequence of length 4, '
                         'containing the values of the corners of the array.')

    initial_condition = create_initial_condition(size, corner_values)

    random_array = iterate_diamond_square(initial_condition, initial_roughness,
                                          roughness_factor, random_seed)

    return random_array
