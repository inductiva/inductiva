"""Utils to create the boundary conditions."""

from abc import ABC, abstractmethod
from typing import Optional, Literal


class BoundaryCondition(ABC):
    """Abstract base class for boundary condition.

    Attributes:
        boundary_name (float): The boundary name of the plate where the boundary
          condition will be apply. 
          The available options are: left, top, right and bottom.
    """

    def __init__(
            self, boundary_name: Literal["left", "top", "right",
                                         "bottom"]) -> None:
        """Initializes a BoundaryCondition object."""
        self.boundary_name = boundary_name

    @abstractmethod
    def to_dict(self):
        """Abstract method to convert the bc properties to a dictionary."""
        pass


class DirichletBC(BoundaryCondition):
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


class NeumannBC(BoundaryCondition):
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
