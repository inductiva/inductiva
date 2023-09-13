"""Utils to check if optional dependencies are installed.

Example usage, e.g., to check if numpy is installed:

extra_msg = "You can install it with `pip install numpy`."

needs_numpy = functools.partial(
    _needs_optional_deps,
    ["numpy"],
    extra_msg,
)

# Now use `needs_numpy` as a decorator for functions that need numpy.
@optional_deps.needs_numpy
def load_from_np():
    ...
"""

import importlib
import functools


def _needs_optional_deps(module_names, extra_msg, func):
    """Decorator to check if optional dependencies are installed."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        missing_imports = []

        for module_name in module_names:
            try:
                importlib.import_module(module_name)
            except ImportError:
                missing_imports.append(module_name)

        if missing_imports:
            missing_imports_quoted = [f"'{m}'" for m in missing_imports]
            missing_imports_str = ", ".join(missing_imports_quoted)

            msg = (f"Dependencies required to use function "
                   f"'{func.__qualname__}' missing: {missing_imports_str}.")

            if extra_msg:
                msg += f" {extra_msg}"

            raise RuntimeError(msg)

        return func(*args, **kwargs)

    return wrapper


def _missing_deps_msg(extra_deps_name):
    return (f"You can install the missing dependencies for '{extra_deps_name}' "
            f"with `pip install 'inductiva[{extra_deps_name}_extra]'`.")


# Apply the two first arguments to _needs_optional_deps function, creating
# a new decorator function with those arguments already set.
needs_molecules_extra_deps = functools.partial(
    _needs_optional_deps,
    ["MDAnalysis", "nglview", "imageio", "matplotlib"],
    _missing_deps_msg("molecules"),
)

needs_fluids_extra_deps = functools.partial(
    _needs_optional_deps,
    ["pyvista", "vtk", "scipy", "xarray", "dask", "imageio", "matplotlib"],
    _missing_deps_msg("fluids"),
)

needs_coastal_extra_deps = functools.partial(
    _needs_optional_deps,
    ["utm", "imageio", "scipy", "matplotlib"],
    _missing_deps_msg("coastal"),
)

needs_common_extra_deps = functools.partial(
    _needs_optional_deps,
    ["imageio", "matplotlib"],
    _missing_deps_msg("common"),
)
