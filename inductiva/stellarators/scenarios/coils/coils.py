"""Class for the coils creation."""
import math
import os
import random
import shutil
import typing
import uuid

from absl import logging
from functools import singledispatchmethod

import numpy as np

from inductiva import scenarios, simulation, stellarators, tasks, types, utils

SIMSOPT_COIL_COEFFICIENTS_FILENAME = 'coil_coefficients.npz'
SIMSOPT_COIL_CURRENTS_FILENAME = 'coil_currents.npz'
SIMSOPT_PLASMA_SURFACE_FILENAME = 'input.final'
SIMSOPT_TEMPLATE_DIR = os.path.join(utils.templates.TEMPLATES_PATH,
                                    'stellarator_coils')
PLASMA_SURFACE_TEMPLATE_FILE_NAME = 'input.example'


class StellaratorCoils(scenarios.Scenario):
    """Represents stellarator coils.

    A stellarator is a magnetic fusion confinement device that possesses a 
    set of external electromagnetic coils that have current running through 
    them. These coils create a magnetic field that is able to confine a plasma, 
    at extremely high temperatures, where the nuclear fusion reactions take 
    place.

    Therefore, designing a stellarator involves defining a set of 
    electromagnetic coils. This set of coils is defined by a small set of 
    NC independent coils, due to the symmetries involved. There is stellarator 
    symmetry, that mirrors this independent set of coils and then 'number of
    field periods (nfp)' rotational symmetry, that repeats this block of coils
    throughout the stellarator in the toroidal direction (the long way around
    the torus). The number of field periods is the number of complete magnetic 
    field repetitions within a device. It represents how many times the magnetic
    field repeats itself along the toroidal direction, lowering the number
    of independent coils that one has to design. Therefore, the total number
    of coils in a stellarator device comes out to 2*NC*nfp.

    For more information about stellarators, there are many articles with 
    very good information. Per Helander's 2014 "Theory of plasma confinement in 
    non-axisymmetric magnetic fields" is a good example. See more details at:
    https://iopscience.iop.org/article/10.1088/0034-4885/77/8/087001

    The class can be initialized normally (directly from __init__) providing
    a list of `Coil` objects or from the class methods so that the individual 
    coils can also be created before creating the StellaratorCoils object. 

    Attributes:
        coils (list): List of `Coil` objects.
        num_field_periods (int): Number of magnetic field periods.
          Refers to the number of complete magnetic field repetitions 
          within a stellarator. Represents how many times the magnetic
          field pattern repeats itself along the toroidal direction.
          Typically varies between 1 and 3.
    """

    valid_simulators = [stellarators.simulators.Simsopt]

    def __init__(self, coils, num_field_periods):
        """Initialize the StellaratorCoils object."""

        self.coils = coils
        self.num_field_periods = num_field_periods

    @classmethod
    def from_circular_curves(cls, num_field_periods, num_coils, coil_currents,
                             major_radius, minor_radius):
        """Create simple circular and equally spaced curves.

        The non-zero coefficients are c0 and c1 for both the Fx and Fy
        Fourier Series and s1 for Fz. This is what makes the coils circular. 

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
                         coil_currents,
                         curves_file,
                         delimiter=','):
        """Create StellaratorCoils from Fourier coefficients loaded from a file.

        This function loads a file containing Fourier coefficients for several 
        coils. The file is expected to have `6*num_coils` many columns, and 
        `order+1` many rows. The columns are in the following order,

        sj_x_coil1,cj_x_coil1,sj_y_coil1,...,sj_x_coil2,cj_x_coil2,...

        Args:
            num_field_periods (int): Number of magnetic field periods.
            coil_currents (list): List of coil currents.
            curves_file (str): Name of the file containing Fourier coefficients.
            delimiter (str): Delimiter used in the file. 

        Returns:
            StellaratorCoils: The created StellaratorCoils instance.
        """

        # Reads all the coefficients from the file to a numpy 2D array
        # with `order+1` columns and `6*num_coils` rows.
        curves_data = np.loadtxt(fname=curves_file,
                                 delimiter=delimiter,
                                 unpack=True)

        # Gets the number of coils.
        num_coils = int((curves_data.shape[0]) / 6)

        # Gets the coefficients for each coil.
        curves_coefficients = np.split(curves_data, num_coils, axis=0)

        coils = [
            Coil(curve_coefficients,
                 coil_current) for curve_coefficients, coil_current in zip(
                     curves_coefficients, coil_currents)
        ]

        return cls(coils, num_field_periods)

    def simulate(
        self,
        simulator: simulation.Simulator = stellarators.simulators.Simsopt(),
        resource_pool_id: typing.Optional[uuid.UUID] = None,
        run_async: bool = False,
        plasma_surface_filepath: typing.Optional[types.Path] = None,
    ) -> tasks.Task:
        """Simulates the scenario.

        Args:
            simulator: The simulator to use for the simulation.
            resource_pool_id: The resource pool to use for the simulation.
            run_async: Whether to run the simulation asynchronously.
            plasma_surface_filepath: Path to the file with the description of
              the plasma surface on which the magnetic field will be calculated.
        """

        if plasma_surface_filepath:
            self.plasma_surface_filepath = utils.files.resolve_path(
                plasma_surface_filepath)
        else:
            logging.info('Plasma surface description not provided. '
                         'Using default file.')
            self.plasma_surface_filepath = os.path.join(
                SIMSOPT_TEMPLATE_DIR, PLASMA_SURFACE_TEMPLATE_FILE_NAME)

        task = super().simulate(
            simulator,
            resource_pool_id=resource_pool_id,
            run_async=run_async,
            coil_coefficients_filename=SIMSOPT_COIL_COEFFICIENTS_FILENAME,
            coil_currents_filename=SIMSOPT_COIL_CURRENTS_FILENAME,
            plasma_surface_filename=SIMSOPT_PLASMA_SURFACE_FILENAME,
            num_field_periods=self.num_field_periods)

        return task

    @singledispatchmethod
    def create_input_files(self, simulator: simulation.Simulator):
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
          defining the coil. Shape: (6, order + 1), where order is the 
          maximum order of the Fourier Series representation.
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


@StellaratorCoils.create_input_files.register
def _(self, simulator: stellarators.simulators.Simsopt, input_dir):  # pylint: disable=unused-argument
    """Creates Simsopt simulation input files."""

    coil_coefficients = [coil.curve_coefficients for coil in self.coils]
    coil_currents = np.array([coil.current for coil in self.coils])

    coil_coefficients_filename = os.path.join(
        input_dir, SIMSOPT_COIL_COEFFICIENTS_FILENAME)
    coil_currents_filename = os.path.join(input_dir,
                                          SIMSOPT_COIL_CURRENTS_FILENAME)

    plasma_surface_filename = os.path.join(input_dir,
                                           SIMSOPT_PLASMA_SURFACE_FILENAME)

    np.savez(coil_coefficients_filename, *coil_coefficients)
    np.savez(coil_currents_filename, coil_currents)
    shutil.copy(self.plasma_surface_filepath, plasma_surface_filename)
