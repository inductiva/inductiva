"""Type definitions."""
import os
import pathlib
from typing import Union, List
from typing_extensions import TypeAlias

Path: TypeAlias = Union[os.PathLike, str]

PathOrStr: TypeAlias = Union[os.PathLike, pathlib.Path, str]

ComputationalResources: TypeAlias = Union["resources.MachineGroup",
                                          "resources.ElasticMachineGroup",
                                          "resources.MPICluster"]

Command: TypeAlias = Union[str, "inductiva.commands.Command"]
Commands: TypeAlias = List[Command]
