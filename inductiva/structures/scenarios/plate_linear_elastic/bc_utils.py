"""Utils to create the different types of boundary conditions."""

from abc import ABC, abstractmethod
from typing import Optional, List, Literal

import json


class BC(ABC):
    """Abstract base class for boundary condition.

    Attributes:
        boundary_name (float): The boundary name of the plate where the boundary
          condition will be apply. 
          The available options are: left, top, right and bottom.
    """

    def __init__(
            self, boundary_name: Literal["left", "top", "right",
                                         "bottom"]) -> None:
        """Initializes a BC object."""
        self.boundary_name = boundary_name

    @abstractmethod
    def to_dict(self):
        """Abstract method to convert the bc properties to a dictionary."""
        pass


class DirichletBC(BC):
    """Dirichlet boundary condition.

    Attributes:
        displacement_x (float): The imposed displacement in x-direction. Only
          fill it if want to resprtit the displamente. 
          Fill this attribute if you want to specify the displacement.
        displacement_y (float): The imposed displacement in y-direction. Only
          fill it if want to resprtit the displamente.
          Fill this attribute if you want to specify the displacement.
    """

    def __init__(self,
                 boundary_name: Literal["left", "top", "right", "bottom"],
                 displacement_x: Optional[float] = None,
                 displacement_y: Optional[float] = None) -> None:
        """Initializes a DirichletBC object."""
        super().__init__(boundary_name)
        self.displacement_x = displacement_x
        self.displacement_y = displacement_y

    def to_dict(self) -> dict:
        """Convert hole properties to a dictionary.

        Returns:
            dict: Hole properties.
        """
        return {
            "boundary_name": self.boundary_name,
            "displacement_x": self.displacement_x,
            "displacement_y": self.displacement_y
        }


class NeumannBC(BC):
    """Neumann boundary condition.

    Attributes:
        tension_x (float): The tension value in x-direction.
        tension_y (float): The tension value in y-direction.
    """

    def __init__(self, boundary_name: Literal["left", "top", "right", "bottom"],
                 tension_x: float, tension_y: float) -> None:
        """Initializes a NeumannBC object."""
        super().__init__(boundary_name)
        self.tension_x = tension_x
        self.tension_y = tension_y

    def to_dict(self) -> dict:
        """Convert hole properties to a dictionary.

        Returns:
            dict: Hole properties.
        """
        return {
            "boundary_name": self.boundary_name,
            "tension_x": self.tension_x,
            "tension_y": self.tension_y
        }


class BoundaryConditionsCase:
    """Boundary conditions for a given case.

    The boundary conditions for a given case are characterized by a list of 
      Dirichlet and Tension boundary conditions.

    Attributes:
        bcs (list): The boundary conditions objects.
    """

    def __init__(self, bcs: List[BC]) -> None:
        """Initializes a BoundaryConditionsCase object."""
        self.bcs = bcs

    def write_to_json(self, json_path: str) -> None:
        """Write the boundary conditions to JSON file.

        Args:
            json_path (str): The JSON file path.
        """

        # Divide the bcs_obj_list into two lists based on their types
        dirichlet_bcs = [bc for bc in self.bcs if isinstance(bc, DirichletBC)]
        neumann_bcs = [bc for bc in self.bcs if isinstance(bc, NeumannBC)]

        # Dirichlet dictionary
        if dirichlet_bcs is not None:
            dirichlet_bcs_dict = []
            for dirichlet_bc in dirichlet_bcs:
                dirichlet_bc_dict = dirichlet_bc.to_dict()
                dirichlet_bcs_dict.append(dirichlet_bc_dict)
            dirichlet_dict = {"dirichlet": dirichlet_bcs_dict}

            # Merge dictionaries: boundary conditions dictionary
            bcs_dict = {**dirichlet_dict}

        # Neumann dictionary
        if neumann_bcs is not None:
            neumann_bcs_dict = []
            for neumann_bc in neumann_bcs:
                neumann_bc_dict = neumann_bc.to_dict()
                neumann_bcs_dict.append(neumann_bc_dict)
            neumann_dict = {"neumann": neumann_bcs_dict}

            # Merge dictionaries: boundary conditions dictionary
            bcs_dict = {**bcs_dict, **neumann_dict}

        # Write JSON file
        with open(json_path, "w", encoding="utf-8") as write_file:
            json.dump(bcs_dict, write_file, indent=4)
