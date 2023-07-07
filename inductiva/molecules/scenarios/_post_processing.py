"""Post process Gromacs simulation outputs."""
from typing import Optional

from absl import logging

import os
import MDAnalysis as mda
import nglview as nv

from uuid import UUID

from inductiva.types import Path
from inductiva.molecules.simulators import GROMACS
from inductiva.utils.templates import (TEMPLATES_PATH,
                                       replace_params_in_template)
from inductiva.tasks import fetch_task_output, get_task_info

import json

SCENARIO_TEMPLATE_DIR = os.path.join(TEMPLATES_PATH, "protein_visualization")
GROMACS_TEMPLATE_INPUT_DIR = "gromacs"


class GROMACSSimulationOutput:
    """Post process Gromacs simulation outputs."""

    def __init__(self, task_id: str, sim_output_path: Path = None):
        """Initializes a `SimulationOutput` object.

        Currently only supports asynchronous simulations."""

        if not get_task_info(task_id)["status"]:
            raise ValueError("Simulation not finished.")

        self.sim_output_dir = sim_output_path

        if sim_output_path is None:
            self.sim_output_dir = os.path.join(os.getcwd(),
                                               f"inductiva_output/{task_id}")

        fetch_task_output(task_id, output_dir=sim_output_path, return_type=None)

        self.template_dir = os.path.join(SCENARIO_TEMPLATE_DIR,
                                         GROMACS_TEMPLATE_INPUT_DIR)

    def render(self,
               trajectory: str = "trajectory.xtc",
               system: str = "Protein-H",
               trr_file="solvated_protein.trr",
               tpr_file="solvated_protein.tpr",
               resource_pool_id: Optional[UUID] = None):
        """Create a movie from the simulation outputs."""

        self.trajectory = trajectory

        commands_path = os.path.join(self.template_dir, "commands.json")
        replace_params_in_template(
            self.template_dir, "visualization_command.json.jinja", {
                "trr_file": trr_file,
                "tpr_file": tpr_file,
                "trajectory": self.trajectory,
                "selected_part": system
            }, commands_path)

        commands = self.read_commands_from_file(commands_path)

        logging.info("Rendering trajectory file .")
        GROMACS().run(self.sim_output_dir,
                      commands=commands,
                      output_dir=self.sim_output_dir,
                      resource_pool_id=resource_pool_id)

        logging.info("Trajectory file rendered.")

    def visualize(self, pdb_file=None):
        """Visualize the simulation outputs in a notebook."""

        if pdb_file is None:
            for filename in os.listdir(self.sim_output_dir):
                if filename.endswith(".pdb"):
                    pdb_file = os.path.join(self.sim_output_dir, filename)

        trajectory = os.path.join(self.sim_output_dir, self.trajectory)
        system = mda.Universe(pdb_file, trajectory)

        view = nv.show_mdanalysis(system)
        view.add_ball_and_stick(
            "all")  # Render the molecules as ball-and-stick models
        view.center()  # Center the view
        view.parameters = {
            "backgroundColor": "white"
        }  # Set the background color

        return view

    def read_commands_from_file(self, commands_path: str):
        "Read list of commands from commands.json file"
        with open(commands_path, "r", encoding="utf-8") as f:
            commands = json.load(f)
        return commands
