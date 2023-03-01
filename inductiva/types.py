"""Type definitions."""
import os
from typing import TypeAlias, Union

Path: TypeAlias = Union[os.PathLike[str], str]
