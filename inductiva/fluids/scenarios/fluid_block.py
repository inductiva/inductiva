"""Describes the physical scenarios and runs its simulation via API."""
import tempfile
import os
import math
from typing import List, Literal, Optional, Union
import xml.etree.ElementTree as ET

from absl import logging

import inductiva_sph
from inductiva_sph import sph_core
import inductiva
from inductiva.fluids._output_post_processing import SimulationOutput
from inductiva.types import Path

# Global variables to define a scenario
TANK_DIMENSIONS = [1, 1, 1]
XML_INPUT_FILENAME = "InputCase.xml"
INPUT_XML_PATH = os.path.join(os.path.dirname(__file__), "xml_files",
                              XML_INPUT_FILENAME)


class FluidBlock:
    """Physical scenario of a general fluid block simulation."""

    def __init__(self,
                 density: float,
                 kinematic_viscosity: float,
                 dimensions: List[float],
                 position: Optional[List[float]] = None,
                 inital_velocity: Optional[List[float]] = None):
        """Initializes a `FluidBlock` object.

        Args:
            density: Density of the fluid in kg/m^3.
            kinematic_viscosity: Kinematic viscosity of the fluid,
                in m^2/s.
            dimensions: A list containing fluid column dimensions,
                in meters.
            position: Position of the fluid column in the tank,
                in meters.
            initial_velocity: Initial velocity of the fluid block
                in the [x, y, z] axes, in m/s.
        """

        self.fluid = sph_core.fluids.FluidProperties(
            density=density, kinematic_viscosity=kinematic_viscosity)

        if len(dimensions) != 3:
            raise ValueError("`fluid_dimensions` must have 3 values.")
        self.dimensions = dimensions

        if position is None:
            self.position = [0.0, 0.0, 0.0]
        else:
            self.position = position

        if inital_velocity is None:
            self.initial_velocity = [0.0, 0.0, 0.0]
        else:
            self.initial_velocity = inital_velocity

    def simulate(self,
                 device: Literal["cpu", "gpu"] = "cpu",
                 engine: Literal["DualSPHysics",
                                 "SPlisHSPlasH"] = "SPlisHSPlasH",
                 simulation_time: float = 1.,
                 particle_radius: float = 0.015,
                 output_dir: Optional[Path] = None,
                 engine_parameters: Union[
                     DualSPHysicsParameters,
                     SPlisHSPlasHParameters] = SPlisHSPlasHParameters):
        """Runs SPH simulation of the Fluid Block scenario.
        Args:
            device: Sets the device for a simulation to be run.
            engine: The software platform to be used for the simulation.
              Available options are (the default is DualSPHysics):
              - SPlisHSPlasH
              - DualSPHysics
            particle_radius: Radius of the discretization particles, in meters.
              Used to control particle spacing. Smaller particle radius means a
              finer discretization, hence more particles.
            time_max: Maximum time of simulation, in seconds.
            output_dir: Directory in which the output files will be saved. If
                not specified, the default directory used for API tasks
                (based on an internal ID of the task) will be used.
            engine_parameters: Simulator specifix parameters.
        """
        self.particle_radius = particle_radius
        self.simulation_time = simulation_time
        self.device = device
        self.output_dir = output_dir
        self.engine_parameters = engine_parameters

        # Create a temporary directory to store simulation input files
        self.input_temp_dir = tempfile.TemporaryDirectory()  #pylint: disable=consider-using-with

        if engine.lower() == "splishsplash" and \
            isinstance(engine_parameters, SPlisHSPlasHParameters):
            sim_output_path = self._splishsplash_simulation()
        elif engine.lower() == "dualsphysics" and \
            isinstance(engine_parameters, DualSPHysicsParameters):
            sim_output_path = self._dualsphysics_simulation()
        else:
            raise ValueError(f"Entered `engine` does not exist or it \
                             does not match with `engine_parameters` class")

        # Delete temporary input directory
        self.input_temp_dir.cleanup()

        inductiva_sph.splishsplash.io_utils.convert_vtk_data_dir_to_netcdf(
            data_dir=os.path.join(sim_output_path, "vtk"),
            output_time_step=self.engine_parameters.output_time_step,
            netcdf_data_dir=os.path.join(sim_output_path, "netcdf"))

        return SimulationOutput(sim_output_path)

    def _splishsplash_simulation(self, input_dir):
        """Runs simulation on SPlisHSPlasH via API.

        Args:
            input_dir: Directory where the input file will be stored.
        """
        # Create a dam break scenario
        scenario = self._create_scenario()

        # Create simulation
        simulation = inductiva_sph.splishsplash.SPlisHSPlasHSimulation(
            scenario=scenario,
            time_max=self.simulation_time,
            particle_radius=self.particle_radius,
            simulation_method=self.simulation_time,
            viscosity_method=self.engine_parameters.viscosity_solver,
            boundary_handling_method=self.engine_parameters.
            boundary_handling_method,
            z_sort=self.engine_parameters.z_sort,
            cfl_method=self.engine_parameters.cfl_method,
            output_time_step=self.engine_parameters.output_time_step,
            output_directory=self.input_temp_dir.name)

        # Create input file
        simulation.create_input_file()
        logging.info("Estimated number of particles %d",
                     self.estimate_num_particles())
        logging.info("Estimated number of time steps %s",
                     math.ceil(self.simulation_time / simulation.time_step))
        logging.info(
            "Number of output time steps %s",
            math.ceil(self.simulation_time /
                      self.engine_parameters.output_time_step))

        sim_output_path = inductiva.sph.splishsplash.run_simulation(
            sim_dir=input_dir.name,
            device=self.device,
            output_dir=self.output_dir)

        return sim_output_path

    def _dualsphysics_simulation(self, input_dir):
        """Runs simulation on DualSPHysics via API.

        Args:
            input_dir: Directory where the input file will be stored.
        """
        # Parse XML file of a dam break scenario
        input_file = ET.parse(INPUT_XML_PATH)
        root = input_file.getroot()

        # Set simulation parameters according to user input
        root.find("./execution/parameters/parameter[@key='TimeMax']").set(
            "value", str(self.simulation_time))
        root.find(".//rhop0").set("value", str(self.fluid.density))
        root.find("./execution/parameters/parameter[@key='Visco']").set(
            "value", str(self.fluid.kinematic_viscosity))
        root.find("./execution/parameters/parameter[@key='TimeOut']").set(
            "value", str(self.engine_parameters.output_time_spet))

        self.update_axis_values_in_xml(
            root=root,
            parameter=".//drawbox[boxfill='solid']/size",
            value=self.fluid_dimensions)
        self.update_axis_values_in_xml(
            root=root,
            parameter=".//drawbox[boxfill='solid']/point",
            value=self.fluid_position)

        particle_size = root.find(".//definition")
        particle_size.set("dp", str(self.particle_radius * 2))

        # Create input file
        input_file.write(os.path.join(input_dir.name, XML_INPUT_FILENAME))

        return inductiva.sph.dualsphysics.run_simulation(
            sim_dir=input_dir.name,
            input_filename=XML_INPUT_FILENAME[:-4],
            device=self.device,
            output_dir=self.output_dir)

    def _create_scenario(self):
        # Create fluid column
        fluid_block = sph_core.fluids.BoxFluidBlock(
            fluid_properties=self.fluid,
            position=self.fluid_position,
            dimensions=self.fluid_dimensions,
            initial_velocity=self.engine_parameters.column_velocity)

        return sph_core.scenarios.DamBreakSPHScenario(
            dimensions=TANK_DIMENSIONS, fluid_blocks=[fluid_block])

    def estimate_num_particles(self):
        """Estimate of the number of SPH particles contained in fluid blocks."""

        # Calculate number of particles for a fluid block
        n_particles_x = round(self.fluid_dimensions[0] /
                              (2 * self.particle_radius)) - 1
        n_particles_y = round(self.fluid_dimensions[1] /
                              (2 * self.particle_radius)) - 1
        n_particles_z = round(self.fluid_dimensions[2] /
                              (2 * self.particle_radius)) - 1

        # Add number of particles to the total sum
        return n_particles_x * n_particles_y * n_particles_z

    def update_axis_values_in_xml(self, root: ET.Element, parameter: str,
                                  value: list):
        """Sets X, Y, Z values of a given parameter in an ElementTree."""
        param = root.find(parameter)
        param.set("x", str(value[0]))
        param.set("y", str(value[1]))
        param.set("z", str(value[2]))
