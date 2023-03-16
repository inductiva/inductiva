"""Classes that define a fluid tank scenario and simulate it via API."""

from absl import logging

from dataclasses import dataclass, field
import math
import tempfile
from typing import List, Literal, Optional

from inductiva_sph import sph_core
from inductiva_sph.splishsplash import SPlisHSPlasHSimulation

from inductiva.types import Path
from inductiva.fluids.shapes import BaseShape
from inductiva.fluids.shapes import Rectangle
from inductiva.fluids.shapes import Circle
from inductiva.fluids.shapes import Cube
from inductiva.fluids.shapes import Cylinder
from inductiva.fluids._fluid_types import WATER
from inductiva.sph import splishsplash


# Tank inlets.
@dataclass
class BaseTankInlet:
    """Base tank inlet."""
    fluid_velocity: float = 1
    position: List[float] = field(default_factory=lambda: [0, 0])


@dataclass
class RectangularTankInlet(BaseTankInlet, Rectangle):
    """Rectangular tank inlet."""
    pass


@dataclass
class CircularTankInlet(BaseTankInlet, Circle):
    """Circular tank inlet."""
    pass


# Tank outlets.
@dataclass
class BaseTankOutlet:
    """Base tank outlet."""
    top_base_position: List[float] = field(default_factory=lambda: [0, 0])


@dataclass
class CubicTankOutlet(BaseTankOutlet, Cube):
    """Cubic tank outlet."""
    pass


@dataclass
class CylindricalTankOutlet(BaseTankOutlet, Cylinder):
    """Cylindrical tank outlet."""
    pass


class FluidTank:
    """Fluid tank."""

    def __init__(
        self,
        shape: BaseShape = Cube(dimensions=[1, 1, 1]),
        fluid: sph_core.fluids.FluidProperties = WATER,
        fluid_level: float = 0,
        inlet: Optional[BaseTankInlet] = CircularTankInlet(radius=0.1,
                                                           position=[0, 0]),
        outlet: Optional[BaseTankOutlet] = CylindricalTankOutlet(
            radius=0.1, height=0.1, top_base_position=[0, 0, -0.1]),
    ):
        self.shape = shape
        self.fluid = fluid
        self.fluid_level = fluid_level
        self.inlet = inlet
        self.outlet = outlet

    def simulate(
        self,
        engine: Literal["SPlisHSPlasH"] = "SPlisHSPlasH",
        device: Literal["CPU", "GPU"] = "CPU",
        output_dir: Optional[Path] = None,
    ):
        """TODO."""

        # Create a temporary directory to store simulation input files
        input_temp_dir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with

        if engine == "SPlisHSPlasH":
            sim_output_path = self._splishsplash_simulation(
                device=device,
                input_dir=input_temp_dir,
                output_dir=output_dir,
            )
        else:
            raise ValueError(f"Invalid engine `{engine}`.")

        # Delete temporary input directory
        input_temp_dir.cleanup()

        # inductiva_sph.splishsplash.io_utils.convert_vtk_data_dir_to_netcdf(
        #     data_dir=os.path.join(sim_output_path, "vtk"),
        #     output_time_step=OUTPUT_TIME_STEP,
        #     netcdf_data_dir=os.path.join(sim_output_path, "netcdf"))

        # return SimulationOutput(sim_output_path)
        return sim_output_path

    def _splishsplash_simulation(self, device, input_dir, output_dir):
        """TODO."""
        # Create a dam break scenario
        scenario = self._create_scenario()

        # Create simulation
        simulation = SPlisHSPlasHSimulation(
            scenario=scenario,
            time_max=5,
            # particle_radius=self.particle_radius,
            # cfl_method=self.cfl_method,
            # output_time_step=OUTPUT_TIME_STEP,
            # viscosity_method=VISCOSITY_SOLVER,
            output_directory=input_dir.name,
        )

        # Create input file
        # simulation.create_input_file()
        # logging.info("Estimated number of particles %d",
        #              self.estimate_num_particles())
        # logging.info("Estimated number of time steps %s",
        #              math.ceil(self.simulation_time / simulation.time_step))
        # logging.info("Number of output time steps %s",
        #              math.ceil(self.simulation_time / OUTPUT_TIME_STEP))

        sim_output_path = splishsplash.run_simulation(sim_dir=input_dir.name,
                                                      device=device,
                                                      output_dir=output_dir)

        return sim_output_path

    # def _create_scenario(self):
    #     scenario = None

    #     # Create fluid column
    #     # fluid_block = sph_core.fluids.BoxFluidBlock(
    #     #     fluid_properties=self.fluid,
    #     #     position=self.fluid_position,
    #     #     dimensions=self.fluid_dimensions,
    #     #     initial_velocity=COLUMN_VELOCITY)

    #     # # Set up scenario
    #     # scenario = sph_core.scenarios.DamBreakSPHScenario(
    #     #     dimensions=TANK_DIMENSIONS, fluid_blocks=[fluid_block])

    #     return scenario