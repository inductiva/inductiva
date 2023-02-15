"""Low-level method that interact with the API for SPH computation.

These functions are to be used inside the higher level constructs of the
Inductiva client. For instance, `scenario.simulate()` uses the `run_simulation`
function.
"""
from .api import invoke_api
from .types import DirPath

# pylint: disable=unused-argument


def run_simulation(sim_dir: DirPath) -> None:
    """Run SplishSplash in the API.

    Args:
        sim_dir: Path to the directory containing the simulation inputs.

    Returns:
        TODO: once we have a class for holding the outputs of the simulation,
        it should be used here
    """
    return invoke_api(locals(), run_simulation)


# pylint: enable=unused-argument
