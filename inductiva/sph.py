"""Low-level method that interact with the API for SPH computation.

These functions are to be used inside the higher level constructs of the
Inductiva client. For instance, `scenario.simulate()` uses the `run_simulation`
function.
"""
from .api import invoke_api

# pylint: disable=unused-argument


def run_simulation(config: dict) -> None:
    """Run generic SPH simulation in the API.

    Args:
        config: JSON object for configuration of the simulation.
            Note: Probably we should create a Type specific for params that
            should be sent serialized as JSON.

    Returns:
        TODO: once we have a class for holding the outputs of the simulation,
        it should be used here
    """
    return invoke_api(locals(), run_simulation)


# pylint: enable=unused-argument
