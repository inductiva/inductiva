"""Molecular Dynamics simulation for water box scenario."""
from functools import singledispatchmethod
from typing import Optional, Literal
import os
import shutil
from uuid import UUID

from inductiva.types import Path
from inductiva.molecules.simulators import GROMACS
from inductiva.simulation import Simulator
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       batch_replace_params_in_template,
                                       replace_params_in_template)
from inductiva.scenarios import Scenario
from inductiva.utils.files import remove_files_with_tag

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "md_water_box")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"


class MDWaterBox(Scenario):
    """Molecular dynamics water box scenario."""

    def __init__(
        self,
        temperature: float = 300,
        box_size: float = 2.3,
    ):
        """The scenario involves simulating a box with water molecules.
        The simulation consists of two main steps: energy minimization
        and molecular dynamics simulation.
        By default, the initial positions of the water molecules are
        arranged uniformally, so to randomize them we perform a
        decorrelation step.
        Args:
            temperature: The temperature of the simulation in Kelvin.
            box_size: The size of the box in nm.
        """
        self.template_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                         GROMACS_TEMPLATE_INPUT_DIR)
        self.temperature = temperature
        if box_size < 2.3:
            raise ValueError("The box size must be greater than 2.3 nm.")
        self.box_size = box_size

    def simulate(
            self,
            simulator: Simulator = GROMACS(),
            output_dir: Optional[Path] = None,
            resource_pool_id: Optional[UUID] = None,
            simulation_time: float = 10,  # ns
            integrator: Literal["md", "sd", "bd"] = "md",
            nsteps_minim: int = 5000):
        """Simulate the water box scenario using molecular dynamics.

        Args:
            output_dir: The output directory to save the simulation results.
            simulation_time: The simulation time in ns.
            integrator: The integrator to use for the simulation. Options:
                - "md" (Molecular Dynamics): Accurate leap-frog algorithm for
                integrating Newton's equations of motion.
                - "sd" (Steepest Descent): Stochastic dynamics integrator with
                leap-frog scheme.
                - "bd" (Brownian Dynamics): Euler integrator for Brownian or
                position Langevin dynamics.

            For more details on the integrators, refer to the GROMACS
            documentation at
            https://manual.gromacs.org/current/user-guide/mdp-options.html.

            nsteps_minim: Number of steps for energy minimization.
        """
        self.nsteps = int(
            simulation_time * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.integrator = integrator
        self.nsteps_minim = nsteps_minim
        commands_path = os.path.join(self.template_dir, "commands.json")
        replace_params_in_template(self.template_dir, "commands.json.jinja",
                                   {"box_size": self.box_size}, commands_path)
        commands = self.read_commands_from_file(commands_path)
        return super().simulate(simulator,
                                output_dir,
                                resource_pool_id=resource_pool_id,
                                commands=commands)

    def simulate_async(
            self,
            simulator: Simulator = GROMACS(),
            resource_pool_id: Optional[UUID] = None,
            simulation_time: float = 10,  # ns
            integrator: Literal["md", "sd", "bd"] = "md",
            nsteps_minim: int = 5000):
        """Simulate the water box scenario using molecular dynamics
        asyncronously.

        Args:
            simulation_time: The simulation time in ns.
            integrator: The integrator to use for the simulation. Options:
                - "md" (Molecular Dynamics): Accurate leap-frog algorithm for
                integrating Newton's equations of motion.
                - "sd" (Steepest Descent): Stochastic dynamics integrator with
                leap-frog scheme.
                - "bd" (Brownian Dynamics): Euler integrator for Brownian or
                position Langevin dynamics.

            For more details on the integrators, refer to the GROMACS
            documentation at
            https://manual.gromacs.org/current/user-guide/mdp-options.html.

            nsteps_minim: Number of steps for energy minimization.
        """
        self.nsteps = int(
            simulation_time * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.integrator = integrator
        self.nsteps_minim = nsteps_minim
        commands_path = os.path.join(self.template_dir, "commands.json")
        replace_params_in_template(self.template_dir, "commands.json.jinja",
                                   {"box_size": self.box_size}, commands_path)
        commands = self.read_commands_from_file(commands_path)
        return super().simulate_async(simulator,
                                      resource_pool_id=resource_pool_id,
                                      commands=commands)

    @singledispatchmethod
    def gen_config(self, simulator: Simulator):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )

    @singledispatchmethod
    def gen_aux_files(self, simulator: Simulator, input_dir: str):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )

    @singledispatchmethod
    def get_config_filename(self, simulator: Simulator):  # pylint: disable=unused-argument
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )


@MDWaterBox.get_config_filename.register
def _(self, simulator: GROMACS):  # pylint: disable=unused-argument
    pass


@MDWaterBox.gen_aux_files.register
def _(self, simulator: GROMACS, input_dir):  # pylint: disable=unused-argument
    """Setup the working directory for the simulation."""
    shutil.copytree(os.path.join(self.template_dir),
                    os.path.join(input_dir),
                    dirs_exist_ok=True)
    remove_files_with_tag(input_dir, ".jinja")


@MDWaterBox.gen_config.register
def _(self, simulator: GROMACS, input_dir):  # pylint: disable=unused-argument
    """Generate the mdp configuration files for the simulation."""
    batch_replace_params_in_template(
        templates_dir=self.template_dir,
        template_filename_paths=[
            "simulation.mdp.jinja",
            "energy_minimization.mdp.jinja",
        ],
        params={
            "integrator": self.integrator,
            "nsteps": self.nsteps,
            "ref_temp": self.temperature,
            "nsteps_minim": self.nsteps_minim,
        },
        output_filename_paths=[
            os.path.join(input_dir, "simulation.mdp"),
            os.path.join(input_dir, "energy_minimization.mdp"),
        ])
