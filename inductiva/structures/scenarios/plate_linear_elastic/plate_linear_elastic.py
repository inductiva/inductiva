# pylint: disable=unused-argument
"""Deformable plate scenario."""

import os
from typing import List, Optional

from functools import singledispatchmethod
from uuid import UUID

from inductiva import structures, simulation, scenarios, tasks, types

import bc_utils
import geometry_utils
import material_utils
import mesh_utils

GEOMETRY_FILENAME = "geometry.json"
MESH_FILENAME = "mesh.msh"
BCS_FILENAME = "bcs.json"
MATERIAL_FILENAME = "material.json"


class DeformablePlate(scenarios.Scenario):
    """Plate linear elastic scenario.

    The plate linear elastic scenario is characterized by the plate, holes,
    boundary conditions, and material properties. The plate and holes together
    form the geometry, which is used to create the mesh for the simulation.
    """

    valid_simulators = [structures.simulators.FEniCSx()]

    def __init__(
        self,
        plate: geometry_utils.RectangularPlate,
        holes: List[geometry_utils.Hole],
        bcs: bc_utils.BoundaryConditionsCase,
        material: material_utils.IsotropicLinearElasticMaterial,
    ):
        """Initializes the plate linear elastic scenario.

        Args:
            plate (RectangularPlate): The rectangular plate geometry.
            holes (Hole): List of holes in the plate.
            bcs (BoundaryConditionsCase): The boundary conditions applied to the
              plate.
            material (IsotropicLinearElasticMaterial): The material properties
              of the plate.
        """
        self.plate = plate
        self.holes = holes
        self.bcs = bcs
        self.material = material
        self.geometry = geometry_utils.GeometricCase(plate_obj=plate,
                                                     list_holes_objs=holes)

    def simulate(
        self,
        simulator: simulation.Simulator = structures.simulators.FEniCSx(),
        resource_pool_id: Optional[UUID] = None,
        run_async: bool = False,
        simulation_time=300,
        output_time_step=10,
    ) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: The simulator to use for the simulation.
            simulation_time: The simulation time, in seconds.
            resource_pool_id: The resource pool to use for the simulation.
            run_async: Whether to run the simulation asynchronously.
        """
        self.simulation_time = simulation_time

        task = super().simulate(simulator,
                                resource_pool_id=resource_pool_id,
                                run_async=run_async)
        return task

    @singledispatchmethod
    def create_input_files(self, simulator: simulation.Simulator):
        pass


@DeformablePlate.create_input_files.register
def _(self,
      simulator: structures.simulators.FEniCSx(),
      input_dir: types.Path) -> None:
    """Creates FEniCSx simulation input files."""

    # Geometry file
    geometry_path = os.join.path(input_dir, GEOMETRY_FILENAME)
    self.geometry.write_to_json(geometry_path)

    # Mesh file
    mesh = mesh_utils.Mesh(self.geometry)
    mesh_path = os.join.path(input_dir, MESH_FILENAME)
    mesh.write_to_msh(mesh_path)

    # BC file
    bcs_path = os.join.path(input_dir, BCS_FILENAME)
    self.bcs.write_to_json(bcs_path)

    # Material file
    material_path = os.join.path(input_dir, MATERIAL_FILENAME)
    self.material.write_to_json(material_path)
