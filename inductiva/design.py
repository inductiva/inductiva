"""Functions to explore design space of simulations."""
import os
from typing import Any, List, Optional
from absl import logging

from inductiva import types
from inductiva.simulation import Simulator
from inductiva.utils import templates
from inductiva.utils import files


def explore_design_space(simulator: Simulator,
                         input_dir: types.Path,
                         template_filename: str,
                         tag: str,
                         values: List[Any],
                         output_dir: Optional[types.Path] = None,
                         **extra_sim_kwargs):
    """Explore over design parameters with a template file.

    Args:
        simulator: Simulator to use.
        input_dir: Directory where the input files (including the template
            file) are located.
        template_filename: Name of the template file.
        tag: Name of the parameter to explore.
        values: Values of the parameter to explore.
        output_dir: Directory where the output files will be stored.
        extra_sim_kwargs: Extra keyword arguments to pass to the simulator.
    """
    input_dir = files.resolve_path(input_dir)
    file_format = os.path.splitext(template_filename)[1]

    if output_dir is None:
        output_dir = input_dir.with_name(f"{input_dir.name}-output")
    else:
        output_dir = files.resolve_path(output_dir)

    logging.info("Exploring design space for attribute \"%s\".", tag)

    for value in values:
        logging.info("Running simulation for %s=%s.", tag, value)

        input_filename = f"input_file_{value}{file_format}"
        input_file = os.path.join(input_dir, input_filename)

        templates.replace_params_in_template(
            str(input_dir),
            template_filename,
            {tag: value},
            input_file,
        )

        simulator.run(
            input_dir=input_dir,
            sim_config_filename=input_filename,
            output_dir=output_dir.joinpath(f"design_{tag}={value}"),
            **extra_sim_kwargs,
        )
