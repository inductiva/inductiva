"""Functions that use the API to compute eigenvalues.

These functions compute the eigenvalues and eigenvectors of
sparse metrices with the corresponding solvers.
"""
import numpy as np
import scipy
from .api import invoke_api
# pylint: disable=unused-argument, disable=redefined-builtin


def eigensolver(matrix: scipy.sparse,
                num_eigenpairs: int) -> tuple[np.ndarray, np.ndarray]:
    return invoke_api(locals(), eigensolver)
