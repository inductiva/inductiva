"""Functions that use the API to compute eigenvalues.

These functions compute the eigenvalues and eigenvectors of
sparse metrices with the corresponding solvers.
"""
import numpy as np
import scipy
from .api import invoke_api
# pylint: disable=unused-argument


# def eigs(matrix: scipy.sparse,
#          num_eigenpairs: int) -> tuple[np.ndarray, np.ndarray]:
#     return invoke_api(locals(), eigs)

def eigs(m: scipy.sparse) -> tuple[np.ndarray, np.ndarray]:
    return invoke_api(locals(), eigs)


# pylint: enable=unused-argument
