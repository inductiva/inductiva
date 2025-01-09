"""DualSPHysics simulator module of the API."""

from typing import List, Optional

from inductiva import types, tasks, simulators


class DualSPHysics(simulators.Simulator):
    """Class to invoke a generic DualSPHysics simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the DualSPHysics simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "dualsphysics"
        self.container_image = self._get_image_uri()

    def run(
        self,
        input_dir: Optional[str],
        shell_script: str,
        *,
        on: types.ComputationalResources,
        storage_dir: Optional[str] = "",
        resubmit_on_preemption: bool = False,
        remote_assets: Optional[List[str]] = None,
        **kwargs,
    ) -> tasks.Task:
        """Executes a DualSPHysics simulation.

        Args:
            input_dir: Directory with simulation input files.
            shell_script: Path to the shell script to run the simulation.
            on: The computational resource to launch the simulation on.
            storage_dir: Directory for storing results.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.

        Returns:
            tasks.Task: An object representing the simulation task.
        """
        commands = [f"bash {shell_script}"]

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           **kwargs)
