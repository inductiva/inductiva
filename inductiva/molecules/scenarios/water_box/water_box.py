"""Water box scenario."""
from functools import singledispatchmethod
from typing import Optional, Literal
import json
import os
import shutil

from inductiva.types import Path
from inductiva.molecules.simulators import GROMACS
from inductiva.simulation import Simulator
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       batch_replace_params_in_template,
                                       replace_params_in_template)
from inductiva.scenarios import Scenario
from inductiva.utils.files import remove_files_with_tag

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "water_box")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"


class WaterBox(Scenario):
    """Water box scenario."""

    def __init__(
        self,
        temperature: float = 300,
        box_size: float = 3.0,
    ):
        """The scenario involves simulating a box filled with water molecules.
        The simulation consists of three main steps: energy minimization, water 
        molecule position decorrelation, and molecular dynamics simulation. 
        Currently, this scenario exclusively supports the GROMACS simulator.
        By default, for a given box size the water molecules are always 
        initially positioned in the same locations following a uniform
        distribution. 
        To introduce randomness in the initial positions of the water molecules, 
        a decorrelation step is performed. This step involves running a short 
        simulation so that their position changes according to a normal 
        distribution.
        Args:
            temperature: The temperature of the simulation in Kelvin.
            box_size: The size of the box in nm.
        """
        self.template_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                         GROMACS_TEMPLATE_INPUT_DIR)
        self.temperature = temperature
        self.box_size = box_size

    def simulate(
            self,
            simulator: Simulator = GROMACS(),
            output_dir: Optional[Path] = None,
            simulation_time: float = 10,  # ns
            integrator: Literal["md", "sd", "bd"] = "md",
            nsteps_minim: int = 5000):
        """Simulate the water box scenario.

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
        replace_params_in_template(
            self.template_dir, "commands.json.jinja",
            {"box_size": self.box_size},
            os.path.join(self.template_dir, "commands.json"))
        commands = self.read_commands_from_file()
        return super().simulate(simulator, output_dir, commands=commands)

    def read_commands_from_file(self):
        "Read list of commands from commands.json file"
        commands_path = os.path.join(self.template_dir, "commands.json")
        with open(commands_path, "r", encoding="utf-8") as f:
            commands = json.load(f)
        return commands

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


@WaterBox.get_config_filename.register
def _(self, simulator: GROMACS):  # pylint: disable=unused-argument
    pass


@WaterBox.gen_aux_files.register
def _(self, simulator: GROMACS, input_dir):  # pylint: disable=unused-argument
    """Setup the working directory for the simulation."""
    shutil.copytree(os.path.join(self.template_dir),
                    os.path.join(input_dir),
                    dirs_exist_ok=True)
    remove_files_with_tag(input_dir, ".jinja")


@WaterBox.gen_config.register
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
