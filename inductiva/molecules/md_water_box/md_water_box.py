"""Molecular Dynamics simulation for water box scenario."""
from functools import singledispatchmethod
from typing import Optional, Literal
import os
import shutil
import io

from inductiva import tasks, resources, simulators, scenarios, utils

SCENARIO_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                     "md_water_box")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"
COMMANDS_TEMPLATE_FILE_NAME = "commands.json.jinja"


class MDWaterBox(scenarios.Scenario):
    """Molecular dynamics water box scenario."""

    valid_simulators = [simulators.GROMACS]

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
            simulator: simulators.Simulator = simulators.GROMACS(),
            machine_group: Optional[resources.MachineGroup] = None,
            run_async: bool = False,
            simulation_time_ns: float = 1,  # ns
            output_timestep_ps: float = 1,  # ps
            integrator: Literal["md", "sd", "bd"] = "md",
            n_steps_min: int = 5000) -> tasks.Task:
        """Simulate the water box scenario using molecular dynamics.

        Args:
            machine_group: The MachineGroup to use for the simulation.
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
            run_async: Whether to run the simulation asynchronously.
        """
        simulator.override_api_method_prefix("mdwater_box")

        self.nsteps = int(
            simulation_time_ns * 1e6 / 2
        )  # convert to fs and divide by the time step of the simulation (2 fs)
        self.integrator = integrator
        self.n_steps_min = n_steps_min
        self.output_frequency = int(
            output_timestep_ps * 1000 /
            2)  # convert to fs and divide by the time step
        # of the simulation (2 fs)

        commands = self.get_commands()
        task = super().simulate(simulator,
                                machine_group=machine_group,
                                commands=commands,
                                run_async=run_async)

        return task

    def get_commands(self):
        """Returns the commands for the simulation."""

        commands_template_path = os.path.join(SCENARIO_TEMPLATE_DIR,
                                              GROMACS_TEMPLATE_INPUT_DIR,
                                              COMMANDS_TEMPLATE_FILE_NAME)

        inmemory_file = io.StringIO()
        utils.templates.replace_params(
            template_path=commands_template_path,
            params={"box_size": self.box_size},
            output_file=inmemory_file,
        )
        commands = self.read_commands_from_file(inmemory_file)

        return commands

    @singledispatchmethod
    def create_input_files(self, simulator: simulators.Simulator):
        pass


@MDWaterBox.create_input_files.register
def _(self, simulator: simulators.GROMACS, input_dir):  # pylint: disable=unused-argument
    """Creates GROMACS simulation input files."""

    template_files_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                      GROMACS_TEMPLATE_INPUT_DIR)

    shutil.copytree(template_files_dir, input_dir, dirs_exist_ok=True)

    utils.templates.batch_replace_params(
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
            "output_frequency": self.output_frequency,
        },
        output_filename_paths=[
            os.path.join(input_dir, "simulation.mdp"),
            os.path.join(input_dir, "energy_minimization.mdp"),
        ],
        remove_templates=True,
    )
