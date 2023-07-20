"""Protein solvation scenario."""

from functools import singledispatchmethod
from typing import Optional, Literal
import os
import shutil
import tempfile
from uuid import UUID

from inductiva.types import Path
from inductiva.molecules.simulators import GROMACS
from inductiva.simulation import Simulator
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       batch_replace_params_in_template,
                                       replace_params_in_template)
from inductiva.scenarios import Scenario
from inductiva.utils import files

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "protein_solvation")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"
COMMANDS_TEMPLATE_FILE_NAME = "commands.json.jinja"


class ProteinSolvation(Scenario):
    """Solvated protein scenario."""

    valid_simulators = [GROMACS]

    def __init__(self, protein_pdb: str, temperature: float = 300):
        """
        Scenario constructor for protein solvation based on the GROMACS
        simulator.
        The three main steps of this scenario are solvation, energy minimization
        and simulation. The user can control the number of steps used to perform
        the energy minimization step and the duration, temperature and
        integrator used to perform the simulation.
        Args:
            protein_pdb: The path to the protein pdb file.
            temperature: The temperature to use for the simulation.
            charged: Whether the protein is charged or not. If None, the charge
            is computed automatically.
        """
        self.template_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                         GROMACS_TEMPLATE_INPUT_DIR)
        self.protein_pdb = files.resolve_path(protein_pdb)
        self.temperature = temperature

    def simulate(
            self,
            simulator: Simulator = GROMACS(),
            output_dir: Optional[Path] = None,
            resource_pool_id: Optional[UUID] = None,
            simulation_time: float = 10,  # ns
            integrator: Literal["md", "sd", "bd"] = "md",
            nsteps_minim: int = 5000,
            visualized_section: str = "Protein-H"):
        """Simulate the solvation of a protein.

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
            visualized_section: The section of the protein to visualize in the
            simulation.
                Options:
                - "Protein-H": The protein with hydrogens. This is the default.
                - "System": The whole system (Protein + Water).
        """

        self.visualized_section = visualized_section

        self.nsteps = int(
            simulation_time * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.integrator = integrator
        self.nsteps_minim = nsteps_minim
        commands = self.get_commands()
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
            nsteps_minim: int = 5000,
            visualized_section: str = "Protein-H"):
        """Simulate the solvation of a protein scenario asyncronously.

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
            visualized_section: The section of the protein to visualize in the
            simulation.
                Options:
                - "Protein-H": The protein with hydrogens. This is the default.
                - "System": The whole system (Protein + Water).
        """

        self.visualized_section = visualized_section

        commands = self.get_commands()

        self.nsteps = int(
            simulation_time * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.integrator = integrator
        self.nsteps_minim = nsteps_minim

        return super().simulate_async(simulator,
                                      resource_pool_id=resource_pool_id,
                                      commands=commands)

    def get_commands(self):
        """Returns the commands for the simulation."""

        commands_template_path = os.path.join(SCENARIO_TEMPLATE_DIR,
                                              GROMACS_TEMPLATE_INPUT_DIR,
                                              COMMANDS_TEMPLATE_FILE_NAME)

        with tempfile.NamedTemporaryFile() as commands_file:
            replace_params_in_template(
                template_path=commands_template_path,
                params={"visualized_section": self.visualized_section},
                output_file_path=commands_file.name,
            )

            commands = self.read_commands_from_file(commands_file.name)

        return commands

    @singledispatchmethod
    def create_input_files(self, simulator: Simulator):
        pass


@ProteinSolvation.create_input_files.register
def _(self, simulator: GROMACS, input_dir):  # pylint: disable=unused-argument
    """Creates GROMACS simulation input files."""

    # rename the pdb file to comply with the naming in the commands list
    shutil.copy(self.protein_pdb, os.path.join(input_dir, "protein.pdb"))

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
            "nsteps_minim": self.nsteps_minim,
        },
        output_filename_paths=[
            os.path.join(input_dir, "simulation.mdp"),
            os.path.join(input_dir, "energy_minimization.mdp"),
        ],
        remove_templates=True,
    )
