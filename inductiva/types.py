"""Type definitions."""
import os
from typing import Union
from typing_extensions import TypeAlias

Path: TypeAlias = Union[os.PathLike, str]
