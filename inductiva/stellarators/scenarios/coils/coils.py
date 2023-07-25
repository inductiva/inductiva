"""Class for the coils creation."""
from inductiva.scenarios import Scenario
import math
import numpy as np


class StellaratorCoils(Scenario):
    """Represents stellarator coils.

    The class can be initialized normally (directly from __init__) providing
    a list of `Coil` objects or from the class methods so that the individual 
    coils can also be created before creating the StellaratorCoils object. 

    Attributes:
        num_coils (int): The number of coils per field period 
        (independent coils).
        coils (list): List of `Coil` objects.
        num_field_periods (int): Number of magnetic field periods.
        Refers to the number of complete magnetic field repetitions 
        within a stellarator. Represents how many times the magnetic
        field pattern repeats itself along the toroidal direction.
    """

    def __init__(self, coils, num_field_periods):
        """Initialize the StellaratorCoils object."""

        self.coils = coils
        self.num_coils = len(coils)
        self.num_field_periods = num_field_periods

    @classmethod
    def from_circular_curves(cls, num_field_periods, num_coils, coil_currents,
                             major_radius, minor_radius):
        """Create simple circular and equally spaced curves.

        Args:
            num_field_periods (int): Number of magnetic field periods.
            num_coils (int): The number of coils per field period.
            coil_currents (list): List of coil currents.
            major_radius (float): distance from the center of the torus 
            (the central axis) to the outer edge of the plasma region.
            minor_radius (float): Radius of the simple circular curves.

        Returns:
            StellaratorCoils: The created StellaratorCoils instance.
        """

        stellarator_symmetry = True
        angle_factor = (1+int(stellarator_symmetry)) * \
            num_field_periods*num_coils

        coils = []
        for i in range(num_coils):
            angle = (i + 0.5) * (2 * np.pi) / (angle_factor)
            curve_coefficients = np.zeros((6, 2))

            # Set the coefficients
            curve_coefficients[1, 0] = math.cos(angle) * major_radius
            curve_coefficients[1, 1] = math.cos(angle) * minor_radius
            curve_coefficients[3, 0] = math.sin(angle) * major_radius
            curve_coefficients[3, 1] = math.sin(angle) * minor_radius
            curve_coefficients[4, 1] = -minor_radius

            # Create the Coil object
            coil = Coil(curve_coefficients, coil_currents[i])

            coils.append(coil)

        # Create StellaratorCoils object
        stellarator_coils = cls(coils, num_field_periods)
        return stellarator_coils

    def simulate(self):
        pass


class Coil:
    """Represents one coil of a stellarator.

    The curve is represented in 3D cartesian coordinates (x, y, z) as a 
    combination of Fourier series [Fx, Fy, Fz], each of which is an expansion 
    over trigonometric functions: 
        Fi = sum_j (sj * sin(j * theta) + cj * cos(j * theta)).
    The curve provided must then be a numpy array with shape `(6, order+1)`, 
    where `order` (the number of columns - 1) is the maximum order of the series 
    (maximum value of j) and the number of rows is 6, one for each of the 
    sj and cj coefficients of each series Fi.

    Attributes: 
        curve_coefficients (np.ndarray): Array with Fourier coefficients 
        defining the coil.
        current (float): Coil current.
    """

    def __init__(self, curve_coefficients, current):
        """Initialize the Coil object."""

        self.curve_coefficients = curve_coefficients
        self.current = current
