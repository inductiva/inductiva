"""Simsopt module of the API."""
from typing import Optional

from inductiva import simulation, tasks, types, resources


class Simsopt(simulation.Simulator):
    """Invokes a simsopt simulation on the API."""

    @property
    def api_method_name(self) -> str:
        return "stellarators.simsopt.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        plasma_surface_filename: str,
        coil_coefficients_filename: str,
        coil_currents_filename: str,
        num_field_periods: int,
        num_iterations: int,
        num_samples: int,
        sigma_scaling_factor: float,
        machine_group: Optional[resources.MachineGroup] = None,
        run_async: bool = False,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
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
            num_iterations: Number of iterations to run for the searching 
              process.
            num_samples: Number of different stellarator samples generated per 
              iteration from a configuration. The samples are generated using
              normal distribution noise for each coefficient of the Fourier 
              Series that describes the coils.
            sigma_scaling_factor: Scaling factor for the sigma value used in 
              the random generation of noise. This argument makes sure that 
              the noise is generated proportionally to each coefficient. It 
              also determines the range of search for new values of the
              coefficients.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            machine_group=machine_group,
            run_async=run_async,
            coil_coefficients_filename=coil_coefficients_filename,
            coil_currents_filename=coil_currents_filename,
            plasma_surface_filename=plasma_surface_filename,
            num_field_periods=num_field_periods,
            num_iterations=num_iterations,
            num_samples=num_samples,
            sigma_scaling_factor=sigma_scaling_factor)
