"""Utilities to create random bathymetric profiles."""

import random

from absl import logging

import numpy as np


def kernel_avg(kernel_offsets,
               bathymetry_array,
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
        max_x, max_y = bathymetry_array.shape
        all_bounded = (src_idxs[:, 0] < max_x) & (src_idxs[:, 1] < max_y)
        # Filter out to make sure all indices are inside the map.
        src_idxs = src_idxs[all_positive & all_bounded]

    # Fetch elements by list of x and list of y coords and return the mean.
    elements = bathymetry_array[src_idxs[:, 0], src_idxs[:, 1]]
    return np.mean(elements)


def diamond_step(bathymetry_array, inter_center_distance, noise_range):
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
    map_length = bathymetry_array.shape[0]
    kernel_scale = inter_center_distance // 2

    for i in range(kernel_scale, map_length, inter_center_distance):
        for j in range(kernel_scale, map_length, inter_center_distance):
            bathymetry_array[i,
                             j] = kernel_avg(diamond_offsets, bathymetry_array,
                                             i, j, kernel_scale)
            bathymetry_array[i, j] += random.uniform(-noise_range, noise_range)


def square_step(bathymetry_array, inter_center_distance, noise_range):
    """The square step performed at a given inter center distance on a map.
                            +-----------+
                            |   | S |   |
                            |---+---+---|
                            | S | X | S |
                            |---+---+---|
                            |   | S |   |
                            +-----------+
    For each square in the array, set the midpoint to be the average of the four
    corner points plus a random value."""

    square_offsets = [[-1, 0], [0, -1], [1, 0], [0, 1]]
    map_length = bathymetry_array.shape[0]
    kernel_scale = inter_center_distance // 2

    # Square Step, rows
    for i in range(kernel_scale, map_length, inter_center_distance):
        for j in range(0, map_length, inter_center_distance):
            bathymetry_array[i,
                             j] = kernel_avg(square_offsets, bathymetry_array,
                                             i, j, kernel_scale)
            bathymetry_array[i, j] += random.uniform(-noise_range, noise_range)

    # Square Step, cols
    for i in range(0, map_length, inter_center_distance):
        for j in range(kernel_scale, map_length, inter_center_distance):
            bathymetry_array[i,
                             j] = kernel_avg(square_offsets, bathymetry_array,
                                             i, j, kernel_scale)
            bathymetry_array[i, j] += random.uniform(-noise_range, noise_range)


def iterate_diamond_square(initial_condition: np.ndarray,
                           initial_roughness: float = 1,
                           roughness_factor: float = 0.5):
    """Performs the diamond-square iteration."""

    if initial_condition.shape[0] != initial_condition.shape[1]:
        logging.error(
            'The initial array must be a square. The current shape is: %f x %f',
            initial_condition.shape[0], initial_condition.shape[1])

    logging.info(
        'Generating the bathymetry with the Diamond-Square Algorithm...')

    map_length = initial_condition.shape[0]
    inter_center_distance = map_length - 1
    roughness = initial_roughness
    random_bathymetry = initial_condition

    while inter_center_distance > 1:

        diamond_step(random_bathymetry, inter_center_distance, roughness)
        square_step(random_bathymetry, inter_center_distance, roughness)

        # Go to a finer grid, and decrease the amount of noise
        inter_center_distance //= 2
        roughness *= roughness_factor

    return random_bathymetry


def generate_inclined_bed_initial_condition(size, max_depth=10):
    """Creates initial conditions for an inclined bed type of profile.

    This is done by setting the left corners of the map lower than the right
    corners. There is some randomness involved to allow for situations where
    the bed is bit rounder,e.g. with more mass on one of the right corners than
    the other."""

    # Start with a flat base shape, all at zero level.
    inclined_bed = np.zeros(shape=(size, size))

    # Define left corners underwater
    inclined_bed[0, 0] = random.uniform(0, max_depth)
    inclined_bed[0, size - 1] = random.uniform(0, max_depth)

    # define right corners above water
    inclined_bed[size - 1, 0] = random.uniform(-max_depth, 0)
    inclined_bed[size - 1, size - 1] = random.uniform(-max_depth, 0)

    return inclined_bed


def create_random_bathymetry(size: int,
                             initial_roughness: float = 1,
                             roughness_factor: float = 0.5,
                             percentile_above_water: float = 20,
                             max_depth: float = 10) -> np.ndarray:
    """Defines a random bathymetry on a regular grid of (x, y) points.

    The random bathymetry is generated with an implementation of the
    Diamond-Square Terrain Generating Algorithm used to create random
    inclined bed profiles.

    The description of the Diamond-Square algorithms can be found in the
    wikipedia (https://en.wikipedia.org/wiki/Diamond-square_algorithm).

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

    Now, it is time for a square step. Observe (below) that there are 4 points,
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
        size: Number of points used to discretize each edge of the domain.
        roughness_factor: Fractional value in range (0, 1) that determines how
          smooth the bathymetry is. Smaller values lead to smoother
          bathymetries.
        percentile_above_water: Percentage of the spatial domain, in range
          (0, 100), that must be above the sea level (<0).
        scale: Scale of the random variables used to initialize the bathymetry.
          For a scale of 10, the initial condition will have random values in
          the range (0, 10) on the left side and (-10,0) on the right side of
          the bathymetry.
    Returns:
        Array with the random bathymetry.
    """

    initial_condition = generate_inclined_bed_initial_condition(size, max_depth)

    random_bathymetry = iterate_diamond_square(initial_condition,
                                               initial_roughness,
                                               roughness_factor)

    # Shift height to ensure that a certain percentage of the domain is above
    # sea level (depth < 0)
    percentile_under_water = np.percentile(random_bathymetry,
                                           percentile_above_water)
    random_bathymetry -= percentile_under_water

    return random_bathymetry
