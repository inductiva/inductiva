"""Protein solvation scenario."""
import os
import shutil
from typing import Optional, Literal
from inductiva import resources, scenarios, simulators, tasks, types, utils

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "protein_solvation")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"


class ProteinSolvation(scenarios.Scenario):
    """Solvated protein scenario."""

    valid_simulators = [simulators.GROMACS]
    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      GROMACS_TEMPLATE_INPUT_DIR)

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

        self.params["protein_pdb"] = utils.files.resolve_path(protein_pdb)
        self.params["ref_temp"] = temperature

    def simulate(
            self,
            simulator: simulators.Simulator = simulators.GROMACS(),
            machine_group: Optional[resources.MachineGroup] = None,
            storage_dir: Optional[types.Path] = "",
            simulation_time_ns: float = 10,  # ns
            output_timestep_ps: float = 1,  # ps
            integrator: Literal["md", "sd", "bd"] = "md",
            n_steps_min: int = 5000,
            ignore_warnings: bool = False) -> tasks.Task:
        """
        Simulate the solvation of a protein.

        Args:
            simulator: The simulator to use for the simulation.
            machine_group: The machine group o use for the simulation.
            run_async: Whether to run the simulation asynchronously.
            storage_dir: The parent directory for storing simulation
            results.
            simulation_time_ns: The simulation time in ns.
            output_timestep_ps: The output timestep in ps.
            integrator: The integrator to use for the simulation. Options:
                - "md" (Molecular Dynamics): Accurate leap-frog algorithm for
                integrating Newton's equations of motion.
                - "sd" (Steepest Descent): Stochastic dynamics integrator with
                leap-frog scheme.
                - "bd" (Brownian Dynamics): Euler integrator for Brownian or
                position Langevin dynamics.
            n_steps_min: Number of steps for energy minimization.
            ignore_warnings: Whether to ignore warnings during grompp (gromacs
            preprocessor). If set to False, the simulation will fail if there
            are warnings. If True, the simulation will run, but the results
            may have inaccuracies. Use with caution.
            n_steps_min: Number of steps for energy minimization.
        """
        simulator.override_api_method_prefix("protein_solvation")

        self.params["nsteps"] = int(
            simulation_time_ns * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.params["output_frequency"] = int(
            output_timestep_ps * 1000 /
            2)  # convert to fs and divide by the time step
        # of the simulation (2 fs)
        self.params["integrator"] = integrator
        self.params["nsteps_minim"] = n_steps_min
        self.params["max_warn"] = 0
        if ignore_warnings:
            self.params["max_warn"] = -1

        commands = self.get_commands()

        task = super().simulate(simulator,
                                machine_group=machine_group,
                                commands=commands,
                                storage_dir=storage_dir)

        return task

    def add_extra_input_files(self, simulator: simulators.GROMACS,
                              input_dir: types.Path):
        """Add protein pdf file into input dir."""

        protein_path = os.path.join(input_dir, "protein.pdb")
        shutil.copy(self.params["protein_pdb"], protein_path)
