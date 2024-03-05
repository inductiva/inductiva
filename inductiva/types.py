"""Type definitions."""
from typing import Union, List, Optional
from typing_extensions import TypeAlias
import pathlib
import os

Path: TypeAlias = Union[os.PathLike, str]

ComputationalResources: TypeAlias = Union["resources.MachineGroup",
                                          "resources.ElasticMachineGroup",
                                          "resources.MPICluster"]

Command: TypeAlias = Union[str, "inductiva.commands.Command"]
Commands: TypeAlias = List[Command]

PathOrStr: TypeAlias = Union[str, pathlib.Path]
OptionalPathOrStr: TypeAlias = Optional[Union[str, pathlib.Path]]
