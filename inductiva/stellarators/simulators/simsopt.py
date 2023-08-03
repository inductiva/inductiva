"""Simsopt module of the API."""
from typing import Optional
from uuid import UUID

from inductiva import types, tasks
from inductiva.simulation import Simulator


class Simsopt(Simulator):
    """Invokes a simsopt simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "stellarators.simsopt.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        magnetic_field_filename: str,
        coil_coefficients_filename: str,
        coil_currents_filename: str,
        plasma_surface_filename: str,
        num_field_periods: int,
        resource_pool_id: Optional[UUID] = None,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            magnetic_field_filename: Name for the output of the simulation.
            coil_coefficients_filename: Name of the file with the Fourier
              Series coefficients of the coils.
            coil_currents_filename: Name of the file with the current in each
              coils.
            plasma_surface_filename: Name of the file with the description of
              the plasma surface on which the magnetic field will be calculated.
            num_field_periods: Number of magnetic field periods.
              Refers to the number of complete magnetic field repetitions 
              within a stellarator. Represents how many times the magnetic
              field pattern repeats itself along the toroidal direction.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            resource_pool_id=resource_pool_id,
            magnetic_field_filename=magnetic_field_filename,
            coil_coefficients_filename=coil_coefficients_filename,
            coil_currents_filename=coil_currents_filename,
            plasma_surface_filename=plasma_surface_filename,
            num_field_periods=num_field_periods)
