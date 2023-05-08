"""Wind tunnel scenario to run an object over a air flow."""

from functools import singledispatchmethod
import os
import shutil
from typing import Optional, List

from inductiva.types import Path
from inductiva.scenarios import Scenario
from inductiva.simulation import Simulator
from inductiva.fluids.simulators import OpenFOAM
from inductiva.utils.templates import batch_replace_params_in_template
from inductiva.utils.files import remove_files_with_tag

OPENFOAM_TEMPLATE_INPUT_DIR = "wind_tunnel_input_template"


class WindTunnel(Scenario):
    """Physical scenario of a general wind tunnel simulation."""

    def __init__(self,
                 object_path: str,
                 flow_velocity: float = 50,
                 domain: Optional[List[float]] = [[-5, 15], [-4, 4], [0, 8]]):

        self.object_path = object_path
        self.flow_velocity = flow_velocity
        self.domain = domain

    def simulate(self,
                 simulator: Simulator = OpenFOAM(),
                 output_dir: Optional[Path] = None,
                 simulation_time: float = 100,
                 write_interval: float = 50,
                 n_cores: int = 1):
        """Simulates the wind tunnel scenario."""

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
    pass


@WindTunnel.gen_config.register
def _(self, simulator: OpenFOAM, input_dir: str):  # pylint: disable=unused-argument
    """Generates the configuration files for OpenFOAM."""

    # The WindTunnel with OpenFOAM requires changing multiple files
    template_dir = os.path.join(os.path.dirname(__file__),
                                "wind_tunnel_input_template/")

    # Copy all files from the template dir to the input directory
    for directory in os.listdir(template_dir):
        shutil.copytree(os.path.join(template_dir, directory),
                        os.path.join(input_dir, directory))

    # Remove all files that have .jinja in input_dir
    remove_files_with_tag(input_dir, ".jinja")

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
