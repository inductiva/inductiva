"""Class for the coils creation."""
from inductiva.scenarios import Scenario
import math
import random
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

        The non-zero coefficients are c_0 and c_1 for both the Fx and Fy
        Fourier Series and s_1 for Fz. This is what makes the coils circular. 

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

        # To spread the coils equally (in a torus) over an angle that
        # is a fraction of 2*pi. 2*num_field_periods*num_coils will be
        # the final number of coils.
        angle_factor = 2 * num_field_periods * num_coils

        coils = []
        for i in range(num_coils):
            toroidal_angle = (i + 0.5) * (2 * np.pi) / (angle_factor)

            # Set the non-zero coefficients
            curve_coefficients = get_circular_curve_coefficients(
                toroidal_angle, major_radius, minor_radius)

            # Create the Coil object
            coil = Coil(curve_coefficients, coil_currents[i])

            coils.append(coil)

        return cls(coils, num_field_periods)

    @classmethod
    def from_random_curves(cls,
                           num_field_periods,
                           num_coils,
                           coil_currents,
                           max_order=12,
                           major_radius=1.0,
                           minor_radius=0.5):
        """Creates random curves.

        Creates a number of 'num_coils' random curves. The random coils are 
        created by first creating simple equally spaced coils and then randomly 
        varying the higher order coefficients of the Fourier series describing 
        the coils. This ensures that the coils are not created on top of each 
        other and are sufficiently separated but with random and complex curves.

        Args:
            num_field_periods (int): Number of magnetic field periods.
            num_coils (int): The number of coils to be randomly created.
            coil_currents (list): List of coil currents.
            max_order (int): Maximum order of the coefficients.
            major_radius (float): distance from the center of the torus 
            (the central axis) to the outer edge of the plasma region.
            minor_radius (float): Radius of the simple initial curves.

        Returns:
            StellaratorCoils: The created StellaratorCoils instance.
        """

        # To spread the coils equally (in a torus) over an angle that
        # is a fraction of 2*pi. 2*num_field_periods*num_coils will be
        # the final number of coils.
        angle_factor = 2 * num_field_periods * num_coils

        coils = []
        for i in range(num_coils):
            toroidal_angle = (i + 0.5) * (2 * np.pi) / (angle_factor)
            curve_coefficients = np.zeros((6, max_order + 1))

            # Base coefficients to make sure the coils are separated.
            curve_coefficients[:, :2] = get_circular_curve_coefficients(
                toroidal_angle, major_radius, minor_radius)

            # Generate random higher coefficients (decreasing the range
            # of the random values to ensure the curve is well behaved).
            for order in range(2, max_order + 1):
                if order == 2:
                    limits = [-0.1, 0.1]

                if 3 <= order <= 5:
                    limits = [-0.05, 0.05]

                if 5 <= order <= 8:
                    limits = [-0.005, 0.005]

                if order > 8:
                    limits = [-0.001, 0.001]

                curve_coefficients[0, order] = random.uniform(*limits)
                curve_coefficients[1, order] = random.uniform(*limits)
                curve_coefficients[2, order] = random.uniform(*limits)
                curve_coefficients[3, order] = random.uniform(*limits)
                curve_coefficients[4, order] = random.uniform(*limits)
                curve_coefficients[5, order] = random.uniform(*limits)

            # Create the Coil object
            coil = Coil(curve_coefficients, coil_currents[i])

            coils.append(coil)

        return cls(coils, num_field_periods)

    @classmethod
    def from_curves_file(cls,
                         num_field_periods,
                         num_coils,
                         coil_currents,
                         curves_file,
                         delimiter=',',
                         max_order=14):
        """Create StellaratorCoils from Fourier coefficients loaded from a file.

        This function loads a file containing Fourier coefficients for several 
        coils. The file is expected to have `6*num_coils` many columns, and 
        `order+1` many rows. The columns are in the following order,

        sin_x_coil1,cos_x_coil1,sin_y_coil1,...,sin_x_coil2,cos_x_coil2,...

        Args:
            num_field_periods (int): Number of magnetic field periods.
            num_coils (int): The number of coils per field period.
            coil_currents (list): List of coil currents.
            curves_file (str): Name of the file containing Fourier coefficients.
            delimiter (str): Delimiter used in the file. 
            max_order (int): Order of the coefficients. 

        Returns:
            StellaratorCoils: The created StellaratorCoils instance.
        """

        coils = []

        # Read the file and extract the coefficients
        with open(curves_file, 'r', encoding='utf-8') as file:
            for curve in range(num_coils):
                # Create a new array to store the coefficients for this curve
                curve_coefficients = np.zeros((6, max_order + 1))

                for order in range(max_order + 1):
                    # Store all the values of the columns corresponding
                    # to the curve
                    line = file.readline().strip().split(delimiter)
                    curve_coefficients[:, order] = [
                        float(val) for val in line[6 * curve:6 * (curve + 1)]
                    ]

                # Create the Coil object
                coil = Coil(curve_coefficients, coil_currents[curve])
                coils.append(coil)

                # Reset the file pointer to the beginning
                file.seek(0)

        return cls(coils, num_field_periods)

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


def get_circular_curve_coefficients(toroidal_angle, major_radius, minor_radius):
    """Sets the non-zero coefficients of a circular curve.

    Args:
        toroidal_angle (float): Angle that defines the position of the
        coil in the torus.
        major_radius (float): distance from the center of the torus 
        (the central axis) to the outer edge of the plasma region.
        minor_radius (float): Radius of the simple initial curves.

    Returns:
        curve_coefficients (np.ndarray): Array with Fourier coefficients
        defining the curve.
    """
    curve_coefficients = np.zeros((6, 2))

    curve_coefficients[1, 0] = math.cos(toroidal_angle) * major_radius
    curve_coefficients[1, 1] = math.cos(toroidal_angle) * minor_radius
    curve_coefficients[3, 0] = math.sin(toroidal_angle) * major_radius
    curve_coefficients[3, 1] = math.sin(toroidal_angle) * minor_radius
    curve_coefficients[4, 1] = -minor_radius

    return curve_coefficients
