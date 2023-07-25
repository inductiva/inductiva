"""Class for the coils creation."""

class StellaratorCoils:
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

    The curve provided must be a list of lists, [sin_x, cos_x, 
    sin_y, cos_y, sin_z, cos_z] where each list has the Fourier 
    coeffiecients wanted.
     
    Attributes: 
        curve (list): List of the Fourier coefficients defining the coil.
        current (float): Coil current.
    """

    def __init__(self, curve, current):
        """Initialize the Coil object."""

        self.curve = curve
        self.current = current
