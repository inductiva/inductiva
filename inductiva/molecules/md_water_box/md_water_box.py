"""Molecular Dynamics simulation for water box scenario."""
from typing import Optional, Literal
import os

from inductiva import tasks, types, resources, simulators, scenarios, utils

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "md_water_box")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"


class MDWaterBox(scenarios.Scenario):
    """Molecular dynamics water box scenario."""

    valid_simulators = [simulators.GROMACS]
    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      GROMACS_TEMPLATE_INPUT_DIR)

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

        self.params["ref_temp"] = temperature
        if box_size < 2.3:
            raise ValueError("The box size must be greater than 2.3 nm.")
        self.params["box_size"] = box_size

    def simulate(
            self,
            simulator: simulators.Simulator = simulators.GROMACS(),
            machine_group: Optional[resources.MachineGroup] = None,
            storage_dir: Optional[types.Path] = "",
            simulation_time_ns: float = 1,  # ns
            output_timestep_ps: float = 1,  # ps
            integrator: Literal["md", "sd", "bd"] = "md",
            n_steps_min: int = 5000) -> tasks.Task:
        """Simulate the water box scenario using molecular dynamics.

        Args:
            simulator: The Simulator to use for the simulation.
            Gromacs is the only supported for now.
            machine_group: The MachineGroup to use for the simulation.
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

            For more details on the integrators, refer to the GROMACS
            documentation at
            https://manual.gromacs.org/current/user-guide/mdp-options.html.
            n_steps_min: Number of steps for energy minimization.
        """
        simulator.override_api_method_prefix("mdwater_box")

        self.params["nsteps"] = int(
            simulation_time_ns * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.params["integrator"] = integrator
        self.params["n_steps_min"] = n_steps_min
        self.params["output_frequency"] = int(
            output_timestep_ps * 1000 /
            2)  # convert to fs and divide by the time step
        # of the simulation (2 fs)

        commands = self.get_commands()

        task = super().simulate(simulator,
                                machine_group=machine_group,
                                commands=commands,
                                storage_dir=storage_dir)

        return task
