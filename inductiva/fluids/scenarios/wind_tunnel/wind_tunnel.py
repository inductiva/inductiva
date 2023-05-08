"""Wind tunnel scenario to run an object over an air flow."""

from functools import singledispatchmethod
import os
import shutil
from typing import Optional, List

from absl import logging

from inductiva.types import Path
from inductiva.scenarios import Scenario
from inductiva.simulation import Simulator
from inductiva.fluids.simulators import OpenFOAM
from inductiva.utils.templates import batch_replace_params_in_template
from inductiva.utils.files import remove_files_with_tag

OPENFOAM_TEMPLATE_INPUT_DIR = "openfoam_template"


class WindTunnel(Scenario):
    """Physical scenario of a general wind tunnel simulation."""

    def __init__(self,
                 object_path: str,
                 flow_velocity: float = 50,
                 domain: Optional[List[float]] = None):
        """Initializes a `WindTunnel` object.
        
        Args:
            object_path: Path to the object that is inserted in the wind tunnel.
            flow_velocity: Velocity of the air flow in m/s.
            domain: List containing the vertices describing the domain
                of the wind tunnel in the [x, y, z] axes, in meters. This is a
                natural description with the default OpenFOAM simulator.
        """

        self.object_path = object_path
        self.flow_velocity = flow_velocity

        if domain is None:
            logging.info("Using a default domain: [[-5, 15], [-4, 4], [0, 8]].")
            self.domain = [[-5, 15], [-4, 4], [0, 8]]

    def simulate(self,
                 simulator: Simulator = OpenFOAM(),
                 output_dir: Optional[Path] = None,
                 simulation_time: float = 100,
                 write_interval: float = 50,
                 n_cores: int = 1):
        """Simulates the wind tunnel scenario.
        
        Args:
            simulator: Simulator to use for the simulation.
            output_dir: Path to the directory where the simulation output
                is downloaded.
            simulation_time: Simulation time, in seconds.
            write_interval: Interval between simulation outputs, in seconds.
            n_cores: Number of cores to use for the simulation.
            """

        self.simulation_time = simulation_time
        self.write_interval = write_interval
        self.n_cores = n_cores

        output_path = super().simulate(
            simulator,
            output_dir=output_dir,
            n_cores=n_cores,
            openfoam_solver="simpleFoam",
        )

        return output_path

    @singledispatchmethod
    @classmethod
    def get_config_filename(cls, simulator: Simulator):  # pylint: disable=unused-argument
        raise ValueError(
            f"Simulator not supported for `{cls.__name__}` scenario.")

    @singledispatchmethod
    def gen_aux_files(self, simulator: Simulator, input_dir: str):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )

    @singledispatchmethod
    def gen_config(self, simulator: Simulator, input_dir: str):
        raise ValueError(
            f"Simulator not supported for `{self.__class__.__name__}` scenario."
        )


@WindTunnel.get_config_filename.register
def _(self, simulator: OpenFOAM):  # pylint: disable=unused-argument
    pass


@WindTunnel.gen_aux_files.register
def _(self, simulator: OpenFOAM, input_dir):  # pylint: disable=unused-argument
    # The WindTunnel with OpenFOAM requires changing multiple files
    template_dir = os.path.join(os.path.dirname(__file__),
                                OPENFOAM_TEMPLATE_INPUT_DIR)

    # Copy all files from the template dir to the input directory
    for directory in os.listdir(template_dir):
        shutil.copytree(os.path.join(template_dir, directory),
                        os.path.join(input_dir, directory))

    # Remove all files that have .jinja in input_dir
    remove_files_with_tag(input_dir, ".jinja")


@WindTunnel.gen_config.register
def _(self, simulator: OpenFOAM, input_dir: str):  # pylint: disable=unused-argument
    """Generates the configuration files for OpenFOAM."""

    # The WindTunnel with OpenFOAM requires changing multiple files
    template_dir = os.path.join(os.path.dirname(__file__),
                                OPENFOAM_TEMPLATE_INPUT_DIR)

    # Add object path to its respective place
    shutil.copy(self.object_path,
                os.path.join(input_dir, "constant/triSurface", "object.obj"))

    batch_replace_params_in_template(
        templates_dir=template_dir,
        template_filename_paths=[
            "0/include/initialConditions_template.openfoam.jinja",
            "system/controlDict_template.openfoam.jinja",
            "system/blockMeshDict_template.openfoam.jinja",
            "system/decomposeParDict_template.openfoam.jinja"
        ],
        params={
            "flow_velocity": self.flow_velocity,
            "simulation_time": self.simulation_time,
            "write_interval": self.write_interval,
            "n_cores": self.n_cores,
            "domain": self.domain,
        },
        output_filename_paths=[
            os.path.join(input_dir, "0/include/initialConditions"),
            os.path.join(input_dir, "system/controlDict"),
            os.path.join(input_dir, "system/blockMeshDict"),
            os.path.join(input_dir, "system/decomposeParDict")
        ])
