"""Classes for scenario simulator mixins."""

import abc
import os
import tempfile

from typing import Optional, Union

from inductiva.fluids.simulators import SPlisHSPlasH
from inductiva.fluids.simulators import DualSPHysics
from inductiva.fluids.simulators import SPlisHSPlasHParameters
from inductiva.fluids.simulators import DualSPHysicsParameters
from inductiva.fluids.post_processing.splishsplash import convert_vtk_data_dir_to_netcdf
from inductiva.types import Path

# TODO: Move this class/submodule elsewhere, to simplify the import.
from inductiva.fluids._output_post_processing import SimulationOutput


class ScenarioSimulatorMixin(abc.ABC):
    """Base class for scenario simulator mixins."""

    def __init__(self):
        """Initializes a `ScenarioSimulatorMixin` object."""
        pass

    def simulate(
        self,
        simulator_params: Union[
            DualSPHysicsParameters,
            SPlisHSPlasHParameters] = SPlisHSPlasHParameters(),
        output_dir: Optional[Path] = None,
    ):
        """Simulates the scenario."""

        with tempfile.TemporaryDirectory() as input_temp_dir:
            output_path = self._simulate(
                simulator_params,
                input_temp_dir,
                output_dir,
            )

        return SimulationOutput(output_path)


class SPlisHSPlasHMixin(ScenarioSimulatorMixin):
    """SPlisHSPlasH mixin.
    
    Defines the methods required by scenarios to be simulated with SPlisHSPlasH.
    """

    def __init__(self):
        """Initializes a `SPlisHSPlasHMixin` object."""
        pass

    def _simulate(self, simulator_params, input_dir, output_dir=None):

        if isinstance(simulator_params, SPlisHSPlasHParameters):
            return self._simulate_splishsplash(simulator_params, input_dir,
                                               output_dir)
        super()._simulate(simulator_params, input_dir, output_dir)

    def _simulate_splishsplash(self, simulator_params, input_dir, output_dir):
        """Simulates the fluid tank with SPlisHSPlasH."""

        self._create_aux_files_splishsplash(input_dir)
        self._replace_params_in_template_splishsplash(input_dir,
                                                      simulator_params)

        simulator = SPlisHSPlasH(sim_dir=input_dir,
                                 sim_config_filename="splishsplash_config.json")

        output_path = simulator.simulate(output_dir=output_dir)

        # TODO: Replace this by a post-processing function, e.g.
        # `on_simulate_end()`?
        convert_vtk_data_dir_to_netcdf(
            data_dir=os.path.join(output_path, "vtk"),
            output_time_step=simulator_params.output_time_step,
            netcdf_data_dir=os.path.join(output_path, "netcdf"))

        return output_path

    def _create_aux_files_splishsplash(self, input_dir):
        """Creates auxiliary files for SPlisHSPlasH simulation."""
        pass

    @abc.abstractmethod
    def _replace_params_in_template_splishsplash(self, input_dir,
                                                 simulator_params):
        """Replaces parameters in a SPlisHSPlasH configuration template."""
        pass


class DualSPHysicsMixin(ScenarioSimulatorMixin):
    """DualSPHysics mixin."""

    def __init__(self):
        """Initializes a `DualSPHysicsMixin` object."""
        pass

    def _simulate(self, simulator_params, input_dir, output_dir=None):

        if isinstance(simulator_params, DualSPHysicsParameters):
            return self._simulate_dualsphysics(simulator_params, input_dir,
                                               output_dir)
        super()._simulate(simulator_params, input_dir, output_dir)

    def _simulate_dualsphysics(self, simulator_params, input_dir, output_dir):
        """Simulates the fluid tank with SPlisHSPlasH."""

        self._create_aux_files_dualsphysics(input_dir)
        self._replace_params_in_template_dualsphysics(input_dir,
                                                      simulator_params)

        simulator = DualSPHysics(sim_dir=input_dir,
                                 sim_config_filename="dualsphysics_config.xml")

        output_path = simulator.simulate(output_dir=output_dir)

        # TODO: Add default post-processing? Convert to NetCDF?

        return output_path

    def _create_aux_files_dualsphysics(self, input_dir):
        """Creates auxiliary files for DualSPHysics simulation."""
        pass

    @abc.abstractmethod
    def _replace_params_in_template_dualsphysics(self, input_dir,
                                                 simulator_params):
        """Replaces parameters in a DualSPHysics configuration template."""
        pass
