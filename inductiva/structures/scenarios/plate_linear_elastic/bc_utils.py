"""Utils to create the different types of boundary conditions."""

from abc import ABC, abstractmethod
from typing import Optional, List, Literal

import json


class BoundaryCoundition(ABC):
    """Abstract base class for boundary condition.

    Attributes:
        boundary_name (float): The boundary name of the plate where the boundary
          condition will be apply. 
          The available options are: left, top, right and bottom.
    """

    def __init__(
            self, boundary_name: Literal["left", "top", "right",
                                         "bottom"]) -> None:
        """Initializes a BoundaryCoundition object."""
        self.boundary_name = boundary_name

    @abstractmethod
    def to_dict(self):
        """Abstract method to convert the bc properties to a dictionary."""
        pass


class DirichletBC(BoundaryCoundition):
    """Dirichlet boundary condition.

    A Dirichlet boundary condition is a fundamental type of boundary condition
    extensively employed in Finite Element Analysis (FEA) to define precise
    values or limitations along certain edges or surfaces within a computational
    domain.

    Mathematically, the Dirichlet boundary condition is utilized to enforce or
    set the solution (such as displacement) at specific boundary points. This is
    done to accurately model real-world restrictions or established behaviors in
    the context of FEA. 

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


class NeumannBC(BoundaryCoundition):
    """A Neumann boundary condition is a vital aspect of Finite Element Analysis
    (FEA), used to specify the flux or the rate of change of a solution variable
    across certain boundaries within a computational domain.

    In mathematical terms, the Neumann boundary condition is employed to
    describe the gradient, flux, or external influence applied to the solution 
    (such as stress) along specific boundary regions. This type of boundary 
    condition is highly valuable in FEA when precise information is available
    about external forces. 

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

    def __init__(self, bcs: List[BoundaryCoundition]) -> None:
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
            list_of_dirichlet_bc_dicts = []
            for dirichlet_bc in dirichlet_bcs:
                dirichlet_bc_dict = dirichlet_bc.to_dict()
                list_of_dirichlet_bc_dicts.append(dirichlet_bc_dict)
            dirichlet_dict = {"dirichlet": list_of_dirichlet_bc_dicts}

            # Merge dictionaries: boundary conditions dictionary
            bcs_dict = {**dirichlet_dict}

        # Neumann dictionary
        if neumann_bcs is not None:
            list_of_neumann_bcs_dicts = []
            for neumann_bc in neumann_bcs:
                neumann_bc_dict = neumann_bc.to_dict()
                list_of_neumann_bcs_dicts.append(neumann_bc_dict)
            neumann_dict = {"neumann": list_of_neumann_bcs_dicts}

            # Merge dictionaries: boundary conditions dictionary
            bcs_dict = {**bcs_dict, **neumann_dict}

        # Write JSON file
        with open(json_path, "w", encoding="utf-8") as write_file:
            json.dump(bcs_dict, write_file, indent=4)
