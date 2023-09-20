"""Coastal area scenario."""

from functools import singledispatchmethod
import math
import os
import shutil
from typing import Literal, Optional
import numpy as np

from absl import logging

from inductiva import tasks, resources, simulators, scenarios, utils
from inductiva import coastal

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "coastal_area")
SWASH_TEMPLATE_SUBDIR = "swash"
SWASH_CONFIG_TEMPLATE_FILENAME = "input.sws.jinja"
SWASH_CONFIG_FILENAME = "input.sws"
SWASH_BATHYMETRY_FILENAME = "bathymetry.bot"


class CoastalArea(scenarios.Scenario):
    """Coastal area scenario.

    This is a simulation scenario for waves propagating in a coastal area
    represented by an arbitrary bathymetric profile (i.e., the map of depths of
    the sea bottom).

    The scenario is simulated in a 2D box (x points east, y points north) with
    dimensions defined by the bathymetric profile. Waves are injected from one
    of the boundaries of the simulation domain with a given amplitude and
    period. The base water level is also configurable.

    The remaining three boundaries of the simulation domain are absorbing.
    Absorption in these boundaries may not be perfect, so small reflections may
    be observed.

    Schematic representation of the simulation scenario: x points right, y
    points up.
    _________________________________
    |                               |
    |                               |
    |                               |
    |  injected waves ->     beach  |
    |                               |
    |                               |
    |_______________________________|
    Note that waves may be injected from other boundaries, and that the beach
    may be located elsewhere in the simulation domain.

    The scenario can be simulated with SWASH.
    """

    valid_simulators = [simulators.SWASH]

    def __init__(
        self,
        bathymetry: coastal.Bathymetry,
        water_level: float = 0,
        wave_source_location: Literal["N", "S", "E", "W"] = "W",
        wave_amplitude: float = 2,
        wave_period: float = 10,
    ):
        """Initializes a `CoastalArea` object.

        Args:
            bathymetry: The bathymetric profile.
            water_level: The water level, in meters.
            wave_source_location: The location of the wave source. Supported
              locations are: N (north), S (south), E (east), W (west),
              corresponding to the upper, lower, right and left boundaries of
              the simulation domain, respectively.
            wave_amplitude: The amplitude of the wave, in meters.
            wave_period: The period of the wave, in seconds.
        """

        if not bathymetry.is_uniform_grid():
            logging.info("The bathymetry is not defined on a uniform grid. "
                         "Attempting to interpolate it to a uniform grid...")
            bathymetry = bathymetry.to_uniform_grid()

        self.bathymetry = bathymetry
        self.water_level = water_level
        self.wave_source_location = wave_source_location
        self.wave_amplitude = wave_amplitude
        self._check_valid_wave()
        self.wave_period = wave_period

    def simulate(
        self,
        simulator: simulators.Simulator = simulators.SWASH(),
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
        simulation_time: float = 100,
        time_step: float = 0.1,
        output_time_step: float = 1,
    ) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: Simulator to use. Supported simulators are: SWASH.
            machine_group: The machine group to use for the simulation.
            simulation_time: Total simulation time, in seconds.
            time_step: Time step, in seconds.
            output_time_step: Time step for the output, in seconds.
        """
        simulator.override_api_method_prefix("coastal_area")

        self.simulation_time = simulation_time
        self.time_step = time_step
        self.output_time_step = output_time_step

        task = super().simulate(
            simulator,
            machine_group=machine_group,
            run_async=run_async,
            sim_config_filename=SWASH_CONFIG_FILENAME,
        )

        return task

    @singledispatchmethod
    def get_config_filename(self, simulator: simulators.Simulator):
        pass

    @singledispatchmethod
    def create_input_files(self, simulator: simulators.Simulator):
        pass

    def _check_valid_wave(self):
        """Checks that the wave amplitude is valid.

        Raises:
            ValueError: If the wave amplitude is not valid.
        """

        depths_grid = self.bathymetry.depths_grid()

        if self.wave_source_location == "W":
            wave_boundary = depths_grid[0, :]

        elif self.wave_source_location == "E":
            wave_boundary = depths_grid[-1, :]

        elif self.wave_source_location == "N":
            wave_boundary = depths_grid[:, -1]

        elif self.wave_source_location == "S":
            wave_boundary = depths_grid[:, 0]

        else:
            raise ValueError("Invalid wave source location. "
                             "Supported locations are: N, S, E, W")

        #The greater the value of depth, the deeper the water, so
        #we want to check the minimum depth using the min value

        minimum_depth = np.min(wave_boundary)
        if self.wave_amplitude > minimum_depth:
            raise ValueError("Wave amplitude is larger than the "
                             "minimum depth of the bathymetry in the "
                             "boundary where waves are generated. \n"
                             "This can be due to a very high wave "
                             "amplitude or the waves are being generated "
                             "on land")


@CoastalArea.get_config_filename.register
def _(cls, simulator: simulators.SWASH):  # pylint: disable=unused-argument
    """Returns the configuration filename for SWASH."""
    return SWASH_CONFIG_FILENAME


@CoastalArea.create_input_files.register
def _(self, simulator: simulators.SWASH, input_dir):  # pylint: disable=unused-argument
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

    # Compute the bathymetry grid spacing.
    bathymetry_x_num = len(self.bathymetry.x_uniques())
    bathymetry_y_num = len(self.bathymetry.y_uniques())

    bathymetry_x_delta = (self.bathymetry.x_ptp()) / bathymetry_x_num
    bathymetry_y_delta = (self.bathymetry.y_ptp()) / bathymetry_y_num

    # SWASH requires the simulation time to be formatted as HHMMSS.sss.
    simulation_time_hmsms = _convert_time_to_hmsms(self.simulation_time)

    # SWASH uses as amplitude the peak-to-peak amplitude.
    wave_amplitude = 2 * self.wave_amplitude

    # All boundaries except the wave source are absorbing.
    absorbing_boundary_locations = ["N", "S", "E", "W"]
    absorbing_boundary_locations.remove(self.wave_source_location)

    utils.templates.replace_params(
        template_path=config_template_file_path,
        params={
            "bathymetry_filename": SWASH_BATHYMETRY_FILENAME,
            "x_range": self.bathymetry.x_range,
            "y_range": self.bathymetry.y_range,
            "x_num": bathymetry_x_num - 1,
            "y_num": bathymetry_y_num - 1,
            "x_delta": bathymetry_x_delta,
            "y_delta": bathymetry_y_delta,
            "water_level": self.water_level,
            "wave_source_location": self.wave_source_location,
            "wave_amplitude": wave_amplitude,
            "wave_period": self.wave_period,
            "absorbing_boundary_locations": absorbing_boundary_locations,
            "simulation_time_hmsms": simulation_time_hmsms,
            "time_step": self.time_step,
            "output_time_step": self.output_time_step,
        },
        output_file=config_file_path,
        remove_template=True,
    )

    bathymetry_file_path = os.path.join(input_dir, SWASH_BATHYMETRY_FILENAME)
    depths_grid = self.bathymetry.depths.reshape(
        (bathymetry_x_num, bathymetry_y_num))

    self.bathymetry.to_bot_file(bathymetry_file_path, depths_grid)


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
