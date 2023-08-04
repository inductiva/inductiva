# pylint: disable=unused-argument
"""Plate linear elastic scenario."""

import os
from typing import Optional, List
from uuid import UUID

from inductiva import types, tasks
from inductiva import structures
from inductiva import simulation
from inductiva import scenarios

from bc_utils import BoundaryConditionsCase
from geometry_utils import CircularHole
from geometry_utils import EllipticalHole
from geometry_utils import GeometricCase
from geometry_utils import RectangularHole
from geometry_utils import RectangularPlate
from material_utils import IsotropicLinearElasticMaterial
from mesh_utils import Mesh


class PlateLinearElastic(scenarios.Scenario):
    """Plate linear elastic scenario.

    The plate linear elastic scenario is characterized by the plate, holes,
    boundary conditions, and material properties. The plate and holes together
    form the geometry, which is used to create the mesh for the simulation.
    """

    valid_simulators = [structures.simulators.FEniCSx]

    def __init__(
        self,
        plate: RectangularPlate,
        holes: List[CircularHole, RectangularHole or EllipticalHole],
        bcs: BoundaryConditionsCase,
        material: IsotropicLinearElasticMaterial,
    ):
        """Initializes the plate linear elastic scenario.

        Args:
            plate: The rectangular plate geometry (RectangularPlate).
            holes: List of holes in the plate (CircularHole, RectangularHole, 
              or EllipticalHole).
            bcs: The boundary conditions applied to the plate.
            material: The material properties of the plate 
              (IsotropicLinearElasticMaterial).
        """
        self.plate = plate
        self.holes = holes
        self.bcs = bcs
        self.material = material

        self.geometry = GeometricCase(plate_obj=plate, list_holes_objs=holes)

    def simulate(
        self,
        simulator: simulation.Simulator = structures.simulators.FEniCSx,
        resource_pool_id: Optional[UUID] = None,
        run_async: bool = False,
        simulation_time=300,
        output_time_step=10,
    ) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: The simulator to use for the simulation.
            simulation_time: The simulation time, in seconds.
            output_time_step: The time step to save the simulation results, in
              seconds.
            resource_pool_id: The resource pool to use for the simulation.
            run_async: Whether to run the simulation asynchronously.
        """
        self.simulation_time = simulation_time
        self.output_time_step = output_time_step

        commands = self.get_commands()

        task = super().simulate(simulator,
                                resource_pool_id=resource_pool_id,
                                run_async=run_async,
                                commands=commands)
        return task


@PlateLinearElastic.create_input_files.register
def _(self,
      simulator: structures.simulators.FEniCSx,
      input_dir: types.Path,
      geometry_filename: str = "geometry.json",
      mesh_filename: str = "mesh.xdmf",
      bcs_filename: str = "bcs.json",
      material_filemane: str = "material.json") -> None:
    """Creates FEniCSx simulation input files."""

    # Geometry file
    geometry_path = os.join.path(input_dir, geometry_filename)
    self.geometry.write_to_json(geometry_path)

    # Mesh file
    mesh = Mesh(self.geometry)
    mesh_path = os.join.path(input_dir, mesh_filename)
    mesh.write_to_xdmf(mesh_path)

    # BC file
    bcs_path = os.join.path(input_dir, bcs_filename)
    self.bcs.write_to_json(bcs_path)

    # Material file
    material_path = os.join.path(input_dir, material_filemane)
    self.material.write_to_json(material_path)
