"""Protein solvation scenario."""
from functools import singledispatchmethod
from typing import Optional
import os
import shutil
from inductiva.types import Path
from inductiva.molecules.simulators import GROMACS
from inductiva.simulation import Simulator
from inductiva.utils.files import resolve_path, get_timestamped_path
from inductiva.utils.misc import split_camel_case
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       batch_replace_params_in_template)
from inductiva.utils.files import remove_files_with_tag

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "protein_solvation")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"


class ProteinSolvation():
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
            working_dir: Optional[Path] = None,
            simulation_time: float = 10,  # ns
            integrator: str = "md",
            nsteps_minim: int = 5000):
        """Simulate the solvation of a protein.
        Args:
            working_dir: The working directory where the simulation 
            will be executed.
            simulation_time: The simulation time in ns. Default is 10 ns.
            integrator: The integrator to use for the simulation. Default is md.
            nsteps_minim: The number of steps to use for the energy minization. 
            Default is 5000.
        """
        self.nsteps = int(
            simulation_time * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.integrator = integrator
        self.nsteps_minim = nsteps_minim
        self.working_dir = self.setup_working_dir(working_dir, self.protein_pdb)
        self.gen_config(simulator)
        self.run_simulation(simulator)

    def setup_working_dir(self, working_dir, protein_pdb):
        """Setup the working directory for the simulation. 
        For this scenario, the input and output directories are the same."""

        # Create working_dir if it does not exist
        if working_dir is None:
            scenario_name_splitted = split_camel_case(self.__class__.__name__)
            working_dir_prefix = "-".join(scenario_name_splitted).lower()
            working_dir = get_timestamped_path(f"{working_dir_prefix}-output")

        working_dir = resolve_path(working_dir)
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)

        # Copy protein pdb to working_dir
        shutil.copy(protein_pdb,
                    os.path.join(working_dir, os.path.basename(protein_pdb)))
        # Copy template files to working_dir
        shutil.copytree(os.path.join(self.template_dir),
                        os.path.join(working_dir),
                        dirs_exist_ok=True)
        # Remove all files that have .jinja in the working_dir
        remove_files_with_tag(working_dir, ".jinja")
        return working_dir

    @singledispatchmethod
    def gen_config(self, simulator: Simulator):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )

    @singledispatchmethod
    def run_simulation(self, simulator: Simulator):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )


@ProteinSolvation.run_simulation.register
def _(self, simulator: GROMACS):
    """Run the simulation using GROMACS."""
    #Solvation
    simulator.run(input_dir=self.working_dir,
                  output_directory=self.working_dir,
                  method_name="pdb2gmx",
                  f=self.protein_pdb,
                  o="protein.gro",
                  water="tip3p",
                  user_input="6")
    simulator.run(input_dir=self.working_dir,
                  output_directory=self.working_dir,
                  method_name="editconf",
                  f="protein.gro",
                  o="protein_box.gro",
                  c="yes",
                  d="1.0",
                  bt="cubic")
    simulator.run(input_dir=self.working_dir,
                  output_directory=self.working_dir,
                  method_name="solvate",
                  cp="protein_box.gro",
                  o="protein_solv.gro",
                  p="topol.top")
    simulator.run(input_dir=self.working_dir,
                  output_directory=self.working_dir,
                  method_name="grompp",
                  f="ions.mdp",
                  c="protein_solv.gro",
                  p="topol.top",
                  o="ions.tpr")
    simulator.run(input_dir=self.working_dir,
                  output_directory=self.working_dir,
                  method_name="genion",
                  s="ions.tpr",
                  o="protein_solv_ions.gro",
                  p="topol.top",
                  pname="NA",
                  nname="CL",
                  neutral="yes")
    # Energy minimization
    simulator.run(input_dir=self.working_dir,
                  output_directory=self.working_dir,
                  method_name="grompp",
                  f="energy_minimization.mdp",
                  c="protein_solv_ions.gro",
                  p="topol.top",
                  o="em.tpr")
    simulator.run(input_dir=self.working_dir,
                  output_directory=self.working_dir,
                  method_name="mdrun",
                  deffnm="em",
                  v="yes")
    # Simulation
    simulator.run(input_dir=self.working_dir,
                  output_directory=self.working_dir,
                  method_name="grompp",
                  f="simulation.mdp",
                  c="em.gro",
                  r="em.gro",
                  p="topol.top",
                  o="solvated_protein.tpr")
    simulator.run(input_dir=self.working_dir,
                  output_directory=self.working_dir,
                  method_name="mdrun",
                  deffnm="solvated_protein",
                  v="yes")


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
