"""Coastal area scenario."""

from functools import singledispatchmethod
import math
import os
from typing import Optional
import shutil
from uuid import UUID

from inductiva.types import Path
from inductiva.scenarios import Scenario
from inductiva.simulation import Simulator
from inductiva.fluids.simulators import SWASH
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       replace_params_in_template)

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "coastal_area")
SWASH_TEMPLATE_SUBDIR = "swash"
SWASH_CONFIG_TEMPLATE_FILENAME = "input.sws.jinja"
SWASH_CONFIG_FILENAME = "input.sws"


class CoastalArea(Scenario):
    """Coastal area scenario."""

    valid_simulators = [SWASH]

    def __init__(self, wave_amplitude: float = 5, wave_period: float = 10):
        """Initializes a `CoastalArea` object.

        Args:
            wave_amplitude: The amplitude of the wave, in meters.
            wave_period: The period of the wave, in seconds.
        """
        self.wave_amplitude = wave_amplitude
        self.wave_period = wave_period

    def simulate(
        self,
        simulator: Simulator = SWASH(),
        output_dir: Optional[Path] = None,
        resource_pool_id: Optional[UUID] = None,
        simulation_time: float = 100,
        time_step: float = 0.1,
        output_time_step: float = 1,
    ):
        """Simulates the scenario.

        Args:
            simulator: Simulator to use.
            output_dir: Directory to store the simulation output.
            resource_pool_id: Resource pool to use for the simulation.
            simulation_time: Total simulation time, in seconds.
            time_step: Time step, in seconds.
            output_time_step: Time step for the output, in seconds.
        """

        self.simulation_time = simulation_time
        self.time_step = time_step
        self.output_time_step = output_time_step

        output_path = super().simulate(
            simulator,
            sim_config_filename=SWASH_CONFIG_FILENAME,
            output_dir=output_dir,
            resource_pool_id=resource_pool_id,
        )

        return output_path

    @singledispatchmethod
    def get_config_filename(self, simulator: Simulator):
        pass

    @singledispatchmethod
    def create_input_files(self, simulator: Simulator):
        pass


@CoastalArea.get_config_filename.register
def _(cls, simulator: SWASH):  # pylint: disable=unused-argument
    """Returns the configuration filename for SWASH."""
    return SWASH_CONFIG_FILENAME


@CoastalArea.create_input_files.register
def _(self, simulator: SWASH, input_dir):  # pylint: disable=unused-argument
    """Creates SPlisHSPlasH simulation input files."""

    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      SWASH_TEMPLATE_SUBDIR)

    shutil.copytree(template_files_dir,
                    input_dir,
                    dirs_exist_ok=True,
                    symlinks=True)

    config_template_file_path, config_file_path = (os.path.join(
        input_dir, file_name
    ) for file_name in [SWASH_CONFIG_TEMPLATE_FILENAME, SWASH_CONFIG_FILENAME])

    replace_params_in_template(
        template_path=config_template_file_path,
        params={
            "wave_amplitude": self.wave_amplitude,
            "wave_period": self.wave_period,
            "simulation_time_hms": _convert_time_to_hmsms(self.simulation_time),
            "time_step": self.time_step,
            "output_time_step": self.output_time_step,
        },
        output_file_path=config_file_path,
        remove_template=True,
    )


def _convert_time_to_hmsms(time: float) -> str:
    """Converts time in seconds to hours, minutes, seconds and milliseconds.

    The output format is HHMMSS.sss, where HH, MM, SS and sss are the
    zero-padded values of hours, minutes, seconds and milliseconds in the time,
    respectively.
    """

    time_m, time_s = divmod(time, 60.)
    time_h, time_m = divmod(time_m, 60.)
    time_ms, time_s = math.modf(time_s)

    return (f"{int(time_h):02d}"
            f"{int(time_m):02d}"
            f"{int(time_s):02d}."
            f"{int(time_ms*1000):03d}")
