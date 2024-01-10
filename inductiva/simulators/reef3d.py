"""Reef3D simulator module of the API."""

from typing import Optional

from inductiva import simulators, types, tasks


class REEF3D(simulators.Simulator):
    """Class to invoke a generic FDS simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "reef3d.reef3d.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[types.Path] = "",
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            on: The computational resource to launch the simulation on. If None
                the simulation is launched in a machine of the default pool.
            other arguments: See the documentation of the base class.
        """
        return super().run(input_dir, on=on, storage_dir=storage_dir)
