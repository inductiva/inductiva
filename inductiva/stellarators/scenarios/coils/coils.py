"""Class for the coils creation."""
from inductiva.scenarios import Scenario


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

    def simulate(self):
        pass


class Coil:
    """Represents one coil of a stellarator.

    The curve is represented in 3D cartesian coordinates (x, y, z) as a 
    combination of Fourier series [Fx, Fy, Fz], each of which is an expansion 
    over trigonometric functions: 
        Fi = sum_j (sj * sin(j * theta) + cj * cos(j * theta)).
    The curve provided must then me a numpy array with dimensions 6*order, 
    where 'order' (the number of columns) is the maximum order of the series 
    (maximum value of j) and the number of rows is 6, one for each of the 
    sj and cj coefficients of each series Fi.
     
    Attributes: 
        curve_coefficients (list): List of the Fourier coefficients defining
        the coil.
        current (float): Coil current.
    """

    def __init__(self, curve_coefficients, current):
        """Initialize the Coil object."""

        self.curve_coefficients = curve_coefficients
        self.current = current
