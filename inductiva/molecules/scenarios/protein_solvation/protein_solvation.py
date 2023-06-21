"""Protein solvation scenario."""
from functools import singledispatchmethod
from typing import Optional, Literal
import os
import shutil

from inductiva.types import Path
from inductiva.molecules.simulators import GROMACS
from inductiva.simulation import Simulator
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       batch_replace_params_in_template)
from inductiva.scenarios import Scenario
from inductiva.utils.files import remove_files_with_tag

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "protein_solvation")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"


class ProteinSolvation(Scenario):
    """Solvated protein scenario."""

    def __init__(
        self,
        protein_pdb: str,
        temperature: float = 300,
    ):
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
        """
        self.template_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                         GROMACS_TEMPLATE_INPUT_DIR)
        self.protein_pdb = protein_pdb
        self.temperature = temperature

    def simulate(
            self,
            simulator: Simulator = GROMACS(),
            output_dir: Optional[Path] = None,
            simulation_time: float = 10,  # ns
            integrator: Literal["md", "sd", "bd"] = "md",
            nsteps_minim: int = 5000):
        """Simulate the solvation of a protein.
        Args:
            output_dir: The  output directory to save the results of the
            simulation. 
            simulation_time: The simulation time in ns. 
            integrator: The integrator to use for the simulation. 
            For more information about the integrator used in the simulation, 
            consult the GROMACS official documentation at 
            https://manual.gromacs.org/current/user-guide/mdp-options.html. 
            nsteps_minim: The number of steps to use for the energy minization. 
        """
        self.nsteps = int(
            simulation_time * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.integrator = integrator
        self.nsteps_minim = nsteps_minim
        return super().simulate(simulator, output_dir)

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


@ProteinSolvation.gen_aux_files.register
def _(self, simulator: GROMACS, input_dir, protein_pdb):  # pylint: disable=unused-argument
    """Setup the working directory for the simulation."""
    # Copy protein pdb to working_dir
    shutil.copy(protein_pdb, os.path.join(input_dir, "protein.pdb"))
    # Copy template files to working_dir
    shutil.copytree(os.path.join(self.template_dir),
                    os.path.join(input_dir),
                    dirs_exist_ok=True)
    # Remove all files that have .jinja in the working_dir
    remove_files_with_tag(input_dir, ".jinja")


@ProteinSolvation.gen_config.register
def _(self, simulator: GROMACS):  # pylint: disable=unused-argument
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
            os.path.join(self.working_dir, "simulation.mdp"),
            os.path.join(self.working_dir, "energy_minimization.mdp"),
        ])
