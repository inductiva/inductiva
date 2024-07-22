"""Simsopt module of the API."""
from typing import Optional

from inductiva import simulators, tasks, types


class SIMSOPT(simulators.Simulator):
    """Invokes a simsopt simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the SIMOPT simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.api_method_name = "stellarators.simsopt.run_simulation"

    def _get_image_uri(self):
        return None

    def run(
        self,
        input_dir: str,
        plasma_surface_filename: str,
        coil_coefficients_filename: str,
        coil_currents_filename: str,
        num_field_periods: int,
        num_iterations: int,
        num_samples: int,
        sigma_scaling_factor: float,
        objectives_weights_filename: str,
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[str] = "",
        resubmit_on_preemption: bool = False,
        extra_metadata: Optional[dict] = None,
        **kwargs,
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
            objectives_weights_filename: Name of the file with the weights for
              each objective function used in the construction of the total
              objective.
            storage_dir: Directory for storing results.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiates with
                `spot=True`.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            objectives_weights_filename=objectives_weights_filename,
            coil_coefficients_filename=coil_coefficients_filename,
            plasma_surface_filename=plasma_surface_filename,
            coil_currents_filename=coil_currents_filename,
            resubmit_on_preemption=resubmit_on_preemption,
            sigma_scaling_factor=sigma_scaling_factor,
            num_field_periods=num_field_periods,
            num_iterations=num_iterations,
            extra_metadata=extra_metadata,
            num_samples=num_samples,
            storage_dir=storage_dir,
            on=on,
            **kwargs,
        )
