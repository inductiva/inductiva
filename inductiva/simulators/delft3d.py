"""Delft3D module of the API."""
from typing import Optional
from inductiva import types, simulators


class Delft3D(simulators.Simulator):
    """Class to invoke a generic Delft3D simulation on the API."""

    _default_delft3d_device = "cpu"

    def __init__(
        self,
        /,
        version: Optional[str] = None,
        use_dev: bool = False,
    ):
        super().__init__(version, use_dev, self._default_delft3d_device)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "delft3d"

    def run(self,
            input_dir: Optional[str],
            *,
            commands: Optional[types.Commands] = None,
            shell_script: Optional[str] = None,
            on: types.ComputationalResources,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[list[str]] = None,
            project: Optional[str] = None,
            **kwargs):

        if commands is None and shell_script is None:
            raise ValueError("Either commands or shell_script "
                             "must be provided.")
        if commands is not None and shell_script is not None:
            raise ValueError("Only one of commands or shell_script "
                             "can be provided.")

        if shell_script is not None:
            #Check if the shell script exists

            self._input_files_exist(input_dir=input_dir,
                                    remote_assets=remote_assets,
                                    shell_script=shell_script)
            commands = [f"bash {shell_script}"]

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           **kwargs)
