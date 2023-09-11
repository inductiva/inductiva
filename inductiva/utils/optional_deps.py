"""Utils to check if optional dependencies are installed.

Example usage, e.g., to check if numpy is installed:

def _check_numpy():
    import numpy

extra_msg = "You can install it with `pip install numpy`."

needs_numpy = functools.partial(
    _needs_optional_deps,
    _check_numpy,
    extra_msg,
)

# Now use `needs_numpy` as a decorator for functions that need numpy.
@optional_deps.needs_numpy
def load_from_np():
    ...
"""

import functools


def _needs_optional_deps(check_fn, extra_info, func):
    """Decorator to check if optional dependencies are installed."""

    def wrapper(*args, **kwargs):
        try:
            check_fn()
        except ImportError as e:
            msg = (f"{str(e)}, required to use function "
                   f"'{func.__qualname__}'.")

            if extra_info:
                msg += f" {extra_info}"

            raise RuntimeError(msg) from e

        return func(*args, **kwargs)

    return wrapper


def _check_molecules_optional_deps():
    # pylint: disable=import-outside-toplevel
    import MDAnalysis as mda
    import nglview as nv

    # unused
    del mda
    del nv
    # pylint: enable=import-outside-toplevel


molecules_missing_deps_msg = (
    "You can install this and other missing dependencies for 'molecules'"
    "with `pip install 'inductiva[molecules_extra]'`.")

# Apply the two first arguments to _needs_optional_deps function, creating
# a new decorator function with those arguments already set.
needs_molecolules_extra_deps = functools.partial(
    _needs_optional_deps,
    _check_molecules_optional_deps,
    molecules_missing_deps_msg,
)
