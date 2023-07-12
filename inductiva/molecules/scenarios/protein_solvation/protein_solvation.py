"""Protein solvation scenario."""

from functools import singledispatchmethod
from typing import Optional, Literal
import os
import shutil
from absl import logging
import numpy as np
from uuid import UUID

from inductiva.types import Path
from inductiva.molecules.simulators import GROMACS
from inductiva.simulation import Simulator
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       batch_replace_params_in_template,
                                       replace_params_in_template)
from inductiva.scenarios import Scenario
from inductiva.utils.files import remove_files_with_tag

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "protein_solvation")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"


class ProteinSolvation(Scenario):
    """Solvated protein scenario."""

    def __init__(self,
                 protein_pdb: str,
                 temperature: float = 300,
                 charged: bool = None):
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
        self.protein_pdb = protein_pdb
        self.temperature = temperature
        self.charged = charged

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

        #Compute charge if it is not provided
        if self.charged is None:
            self.charged = self.compute_charge(
                resource_pool_id=resource_pool_id)

        #Edit commands.json according to the charge of the protein
        commands_path = os.path.join(self.template_dir, "commands.json")

        replace_params_in_template(
            self.template_dir, "commands.json.jinja", {
                "pdb_file": self.protein_pdb,
                "charged": self.charged,
                "visualized_section": visualized_section
            }, commands_path)

        commands = self.read_commands_from_file(commands_path)
        self.nsteps = int(
            simulation_time * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.integrator = integrator
        self.nsteps_minim = nsteps_minim
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

        #Compute charge if it is not provided
        if self.charged is None:
            self.charged = self.compute_charge(
                resource_pool_id=resource_pool_id)

        #Edit commands.json according to the charge of the protein
        commands_path = os.path.join(self.template_dir, "commands.json")

        replace_params_in_template(
            self.template_dir, "commands.json.jinja", {
                "pdb_file": self.protein_pdb,
                "charged": self.charged,
                "visualized_section": visualized_section
            }, commands_path)

        commands = self.read_commands_from_file(commands_path)

        self.nsteps = int(
            simulation_time * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.integrator = integrator
        self.nsteps_minim = nsteps_minim

        self.task_id = super().simulate_async(simulator,
                                              resource_pool_id=resource_pool_id,
                                              commands=commands)

        return self.task_id

    def compute_charge(self,
                       simulator: Simulator = GROMACS(),
                       resource_pool_id: Optional[UUID] = None):
        """Check if the protein is charged.

        If the protein is not charged, the charge is computed by reading
        the topology file generated by the GROMACS tool gmx pdb2gmx."""

        logging.info("Computing the charge of the protein")

        protein_directory = os.path.dirname(self.protein_pdb)

        commands_path = os.path.join(self.template_dir,
                                     "charge_computation.json")
        protein_pdb = os.path.basename(self.protein_pdb)
        replace_params_in_template(self.template_dir,
                                   "charge_computation.json.jinja",
                                   {"pdb_file": protein_pdb}, commands_path)

        commands = self.read_commands_from_file(commands_path)

        simulator.run(protein_directory,
                      commands=commands,
                      output_dir=protein_directory,
                      resource_pool_id=resource_pool_id)

        topology_file = os.path.join(protein_directory, "topol.top")

        #Information about the charge of each residue is stored in the topology
        #file in the lines starting with "; residue". Summing the charges of
        #all residues gives the total charge of the protein.

        with open(topology_file, "r", encoding="utf-8") as file:
            charge = np.sum([
                float(line.split()[-1])
                for line in file
                if line.startswith("; residue")
            ])

        is_charged = abs(charge) > 1e-6
        return is_charged

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


@ProteinSolvation.get_config_filename.register
def _(self, simulator: GROMACS):  # pylint: disable=unused-argument
    pass


@ProteinSolvation.gen_aux_files.register
def _(self, simulator: GROMACS, input_dir):  # pylint: disable=unused-argument
    """Setup the working directory for the simulation."""
    # rename the pdb file to comply with the naming in the commands list
    shutil.copy(self.protein_pdb, os.path.join(input_dir, "protein.pdb"))
    shutil.copytree(os.path.join(self.template_dir),
                    os.path.join(input_dir),
                    dirs_exist_ok=True)
    remove_files_with_tag(input_dir, ".jinja")


@ProteinSolvation.gen_config.register
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
