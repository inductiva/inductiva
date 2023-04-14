"""Base class for scenarios."""

from abc import ABC
import tempfile
from typing import Optional

from inductiva.types import Path
from inductiva.simulation import Simulator
from inductiva.utils.files import resolve_path, get_timestamped_path
from inductiva.utils.misc import split_camel_case


class Scenario(ABC):
    """Base class for scenarios."""

    def simulate(
        self,
        simulator: Simulator,
        output_dir: Optional[Path] = None,
        **kwargs,
    ):
        """Simulates the scenario."""
        if output_dir is None:
            scenario_name_splitted = split_camel_case(self.__class__.__name__)
            output_dir_prefix = "-".join(scenario_name_splitted).lower()
            output_dir = get_timestamped_path(f"{output_dir_prefix}-output")

        output_dir = resolve_path(output_dir)

        with tempfile.TemporaryDirectory() as input_dir:

            self.gen_aux_files(simulator, input_dir)
            self.gen_config(simulator, input_dir)

            args = ()
            config_filename = self.get_config_filename(simulator)
            if config_filename:
                args += (config_filename,)

            output_path = simulator.run(
                input_dir,
                *args,
                output_dir=output_dir,
                **kwargs,
            )

        return output_path
