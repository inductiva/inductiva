"""Utility functions for the demo scripts."""
import numpy as np
import scipy


def get_square_tridiagonal_h_matrix(size):
    diags = [
        np.random.normal(size=size - 1),
        np.random.normal(size=size),
        np.random.normal(size=size - 1),
    ]

    m = scipy.sparse.diags(diagonals=diags, offsets=[-1, 0, 1], format="csr")

    return m + m.transpose()
