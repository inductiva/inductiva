"""Base class for scenarios."""

from abc import ABC
import tempfile
from typing import Optional, List

from inductiva.types import Path
from inductiva.simulation import Simulator, Command
from inductiva.utils.files import resolve_path, get_timestamped_path
from inductiva.utils.misc import split_camel_case


class Scenario(ABC):
    """Base class for scenarios."""

    def _setup_output_dir(self, output_dir: str):
        """Setup the scenario output directory."""
        if output_dir is None:
            scenario_name_splitted = split_camel_case(self.__class__.__name__)
            output_dir_prefix = "-".join(scenario_name_splitted).lower()
            output_dir = get_timestamped_path(f"{output_dir_prefix}-output")
        output_dir = resolve_path(output_dir)
        return output_dir

    def _setup_config(self, simulator: Simulator, input_dir: Path):
        """Setup the scenario configuration files and arguments."""
        self.gen_aux_files(simulator, input_dir)
        self.gen_config(simulator, input_dir)

        args = ()
        config_filename = self.get_config_filename(simulator)
        if config_filename:
            args += (config_filename,)
        return args

    def gen_pipeline(self, simulator: Simulator):  # pylint: disable=unused-argument
        """Generate the pipeline for the scenario. To be implemented in subclasses."""
        return None

    def simulate(
        self,
        simulator: Simulator,
        output_dir: Optional[Path] = None,
        **kwargs,
    ):
        """Simulates the scenario for a single simulator call."""
        output_dir = self._setup_output_dir(output_dir)

        with tempfile.TemporaryDirectory() as input_dir:
            args, pipeline = self._setup_config(simulator, input_dir)
            pipeline = self.gen_pipeline(simulator)
            if pipeline is None:
                output_path = simulator.run(
                    input_dir,
                    *args,
                    output_dir=output_dir,
                    **kwargs,
                )
            else:
                output_path = simulator.run_pipeline(input_dir,
                                                     *args,
                                                     output_dir=output_dir,
                                                     pipeline=pipeline)
        return output_path

    def simulate_async(
        self,
        simulator: Simulator,
        **kwargs,
    ):
        """Simulates the scenario asychronously."""

        with tempfile.TemporaryDirectory() as input_dir:
            args = self._setup_config(simulator, input_dir)
            task_id = simulator.run_async(
                input_dir,
                *args,
                **kwargs,
            )
        return task_id
