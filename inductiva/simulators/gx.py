"""Gx module of the API for numerical simulations."""
from typing import List, Optional

from inductiva import types, tasks, simulators


class Gx(simulators.Simulator):
    """Class to invoke a generic Gx simulation on the API.
    """

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """Initialize the Gx simulator.

        Args:
            version (str): The version of the simulator to use. If None, the
                latest available version in the platform is used.
            use_dev (bool): Request use of the development version of
                the simulator. By default (False), the production version
                is used.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "gx"
        self.container_image = self._get_image_uri()

    def run(self,
            input_dir: Optional[str],
            sim_config_filename: str,
            *,
            on: types.ComputationalResources,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[List[str]] = None,
            **kwargs) -> tasks.Task:
        """Run the simulation.
        Args:
            input_dir: Path to the directory of the simulation input files.
            on: The computational resource to launch the simulation on.
            sim_config_filename: The name of the simulation configuration file.
            other arguments: See the documentation of the base class.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous execution attempts were preempted. Only applicable when
                using a preemptible resource, i.e., resource instantiated with
                `spot=True`.
            remote_assets: Additional remote files that will be copied to
                the simulation directory.
        """

        commands = [f"gx {sim_config_filename}"]

        return super().run(input_dir,
                           on=on,
                           storage_dir=storage_dir,
                           commands=commands,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           **kwargs)
