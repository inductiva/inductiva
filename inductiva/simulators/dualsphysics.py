"""DualSPHysics simulator module of the API."""

from typing import Optional

from inductiva import types, tasks, simulators


class DualSPHysics(simulators.Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "sph.dualsphysics.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        commands: types.Commands,
        on: Optional[types.ComputationalResources] = None,
        storage_dir: Optional[types.Path] = "",
        extra_metadata: Optional[dict] = None,
    ) -> tasks.Task:
        """Executes a DualSPHysics simulation.

        Args:
            input_dir: Directory with simulation input files.
            sim_config_filename: Simulation config file.
            on: The computational resource to launch the simulation on. If None
                the simulation is submitted to a machine in the default pool.
            storage_dir: Directory for storing results.

        Returns:
            tasks.Task: An object representing the simulation task.
        """
        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           extra_metadata=extra_metadata)
