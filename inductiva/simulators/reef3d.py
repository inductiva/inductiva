"""Reef3D simulator module of the API."""

from typing import Optional

from inductiva import simulators, types, tasks, resources


class REEF3D(simulators.Simulator):
    """Class to invoke a generic FDS simulation on the API."""

    def __init__(self):
        super().__init__()
        self.api_method_name = "reef3d.reef3d.run_simulation"

    def run(
        self,
        input_dir: types.Path,
        machine_group: Optional[resources.MachineGroup] = None,
        storage_dir: Optional[types.Path] = "",
        save_sim_metadata_dir: Optional[types.Path] = None,
        extra_sim_metadata_to_save: Optional[dict] = None,
    ) -> tasks.Task:
        """Run the simulation.

        Args:
            input_dir: Path to the directory of the simulation input files.
            sim_config_filename: Name of the simulation configuration file.
            other arguments: See the documentation of the base class.
        """
        return super().run(
            input_dir,
            machine_group=machine_group,
            storage_dir=storage_dir,
            save_sim_metadata_dir=save_sim_metadata_dir,
            extra_sim_metadata_to_save=extra_sim_metadata_to_save)
