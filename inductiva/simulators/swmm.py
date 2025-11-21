"""SWMM module of the API."""
from typing import Optional, Union

from inductiva import simulators, tasks, types


class SWMM(simulators.Simulator):
    """Class to invoke a generic SWMM simulation on the API."""

    def __init__(self, /, version: Optional[str] = None, use_dev: bool = False):
        """
        Initialize the SWMM simulator wrapper.

        Args:
            version (str): Simulator version to use. If None, the latest
                available version on the platform is used.
            use_dev (bool): Request the development image instead of production.
        """
        super().__init__(version=version, use_dev=use_dev)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "swmm"

    def run(
        self,
        input_dir: Optional[str],
        commands: types.Commands,
        *,
        on: types.ComputationalResources,
        storage_dir: Optional[str] = "",
        resubmit_on_preemption: bool = False,
        remote_assets: Optional[Union[str, list[str]]] = None,
        project: Optional[str] = None,
        time_to_live: Optional[str] = None,
        on_finish_cleanup: Optional[Union[str, list[str]]] = None,
        **kwargs,
    ) -> tasks.Task:
        """Run a list of SWMM commands.

        Args:
            input_dir: Path to the directory containing the input files.
            commands: List of commands to run inside the container
                (for example ["swmm5 model.inp swmm.rpt"]).
            on: The computational resource to launch the simulation on.
            storage_dir: Directory for storing simulation results.
            resubmit_on_preemption (bool): Resubmit task for execution when
                previous attempts were preempted (applicable to spot resources).
            remote_assets: Additional remote files copied into the simulation.
            project: Project name to associate the task with.
            time_to_live: Maximum allowed runtime for the task.
            on_finish_cleanup: Cleanup commands/scripts executed after finish.
            kwargs: Additional keyword arguments passed to the base run method.
        """
        return super().run(
            input_dir,
            on=on,
            commands=commands,
            storage_dir=storage_dir,
            resubmit_on_preemption=resubmit_on_preemption,
            remote_assets=remote_assets,
            project=project,
            time_to_live=time_to_live,
            on_finish_cleanup=on_finish_cleanup,
            **kwargs,
        )
