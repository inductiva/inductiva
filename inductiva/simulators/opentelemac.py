"""OpenTelemac module of the API."""
from typing import Optional, Union

from inductiva import simulators, types


class OpenTelemac(simulators.Simulator):
    """Class to invoke a generic OpenTelemac simulation on the API."""

    _default_opentelemac_device = "cpu"

    def __init__(
        self,
        /,
        version: Optional[str] = None,
        use_dev: bool = False,
    ):
        super().__init__(version, use_dev, self._default_opentelemac_device)
        self.simulator = "arbitrary_commands"
        self.simulator_name_alias = "opentelemac"

    def run(self,
            input_dir: Optional[str],
            commands: types.Commands,
            *,
            on: types.ComputationalResources,
            storage_dir: Optional[str] = "",
            resubmit_on_preemption: bool = False,
            remote_assets: Optional[Union[str, list[str]]] = None,
            project: Optional[str] = None,
            time_to_live: Optional[str] = None,
            **kwargs):

        return super().run(input_dir,
                           on=on,
                           commands=commands,
                           storage_dir=storage_dir,
                           resubmit_on_preemption=resubmit_on_preemption,
                           remote_assets=remote_assets,
                           project=project,
                           time_to_live=time_to_live,
                           **kwargs)
