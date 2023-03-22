"""Low-level method that interacts with the API for SPlisHSPlasH computation.

These functions are to be used inside the higher level constructs of the
Inductiva client. For instance, `scenario.simulate()` uses the `run_simulation`
function.
"""
import pathlib
from typing import Literal, Optional

from inductiva.api import invoke_api
from inductiva.types import Path

# pylint: disable=unused-argument


def run_simulation(sim_dir: Path,
                   input_filename: str = "splishsplash_input.json",
                   device: Literal["gpu", "cpu"] = "cpu",
                   output_dir: Optional[Path] = None) -> pathlib.Path:
    """Run SplishSplash in the API.

    Args:
        sim_dir: Path to the directory containing the simulation inputs.

    Returns:
        TODO(https://github.com/inductiva/inductiva/issues/99): Add class
            for holding the outputs of the simulation
    """
    params = locals()
    del params["output_dir"]

    return invoke_api(params, run_simulation, output_dir=output_dir)


# pylint: enable=unused-argument
