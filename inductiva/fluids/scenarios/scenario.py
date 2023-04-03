"""Base class for scenarios."""

from abc import ABC
import tempfile
from typing import Optional

from inductiva.types import Path
from inductiva.simulation import Simulator


class Scenario(ABC):
    """Base class for scenarios."""

    def simulate(
        self,
        simulator: Simulator,
        output_dir: Optional[Path] = None,
        **kwargs,
    ):
        """Simulates the scenario."""

        simulator_type = type(simulator)

        try:
            sim_config_filename = self._config_filename[simulator_type]
            gen_aux_fn = self._gen_aux[simulator_type]
            gen_config_fn = self._gen_config[simulator_type]
        except KeyError:
            raise NotImplementedError(
                f"Scenario `{type(self).__name__}` does not support "
                f"simulator `{simulator_type.__name__}`.") from None

        with tempfile.TemporaryDirectory() as input_dir:

            gen_aux_fn(self, input_dir)
            gen_config_fn(self, input_dir)

            output_path = output_path = simulator.run(
                input_dir,
                sim_config_filename,
                output_dir=output_dir,
                **kwargs,
            )

        return output_path
