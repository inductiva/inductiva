"""Molecular Dynamics simulation for water box scenario."""
from functools import singledispatchmethod
from typing import Optional, Literal
import os
import shutil
import tempfile

from inductiva import tasks, resources
from inductiva.molecules.simulators import GROMACS
from inductiva.simulation import Simulator
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       batch_replace_params_in_template,
                                       replace_params_in_template)
from inductiva.scenarios import Scenario
from .post_processing import MDWaterBoxOutput

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "md_water_box")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"
COMMANDS_TEMPLATE_FILE_NAME = "commands.json.jinja"


class MDWaterBox(Scenario):
    """Molecular dynamics water box scenario."""

    valid_simulators = [GROMACS]

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

        self.temperature = temperature
        if box_size < 2.3:
            raise ValueError("The box size must be greater than 2.3 nm.")
        self.box_size = box_size

    def simulate(
            self,
            simulator: Simulator = GROMACS(),
            machine_group: Optional[resources.MachineGroup] = None,
            run_async: bool = False,
            simulation_time_ns: float = 10,  # ns
            integrator: Literal["md", "sd", "bd"] = "md",
            n_steps_min: int = 5000) -> tasks.Task:
        """Simulate the water box scenario using molecular dynamics.

        Args:
            machine_group: The MachineGroup to use for the simulation.
            simulation_time_ns: The simulation time in ns.
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

            machine_group: The machine group to use for the simulation.
            n_steps_min: Number of steps for energy minimization.
            run_async: Whether to run the simulation asynchronously.
        """
        self.nsteps = int(
            simulation_time_ns * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.integrator = integrator
        self.n_steps_min = n_steps_min

        commands = self.get_commands()
        task = super().simulate(simulator,
                                machine_group=machine_group,
                                commands=commands,
                                run_async=run_async)

        task.set_output_class(MDWaterBoxOutput)

        return task

    def get_commands(self):
        """Returns the commands for the simulation."""

        commands_template_path = os.path.join(SCENARIO_TEMPLATE_DIR,
                                              GROMACS_TEMPLATE_INPUT_DIR,
                                              COMMANDS_TEMPLATE_FILE_NAME)

        with tempfile.NamedTemporaryFile() as commands_file:
            replace_params_in_template(
                template_path=commands_template_path,
                params={"box_size": self.box_size},
                output_file_path=commands_file.name,
            )

            commands = self.read_commands_from_file(commands_file.name)

        return commands

    @singledispatchmethod
    def create_input_files(self, simulator: Simulator):
        pass


@MDWaterBox.create_input_files.register
def _(self, simulator: GROMACS, input_dir):  # pylint: disable=unused-argument
    """Creates GROMACS simulation input files."""

    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      GROMACS_TEMPLATE_INPUT_DIR)

    shutil.copytree(template_files_dir, input_dir, dirs_exist_ok=True)

    batch_replace_params_in_template(
        templates_dir=input_dir,
        template_filenames=[
            "simulation.mdp.jinja",
            "energy_minimization.mdp.jinja",
        ],
        params={
            "integrator": self.integrator,
            "nsteps": self.nsteps,
            "ref_temp": self.temperature,
            "nsteps_minin": self.n_steps_min,
        },
        output_filename_paths=[
            os.path.join(input_dir, "simulation.mdp"),
            os.path.join(input_dir, "energy_minimization.mdp"),
        ],
        remove_templates=True,
    )
