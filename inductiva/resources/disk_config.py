"""
This module contains the DiskConfig class which is used to store the
configuration for the disk.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass(frozen=True)
class DiskConfig:
    """Dataclass to represent the disk configuration.
    
    A disk has a size in GB and can be resized based on the free space left.
    To enable resizing, the following attributes must be provided:
    resize_trigger_gb, resize_increment_gb, and max_size_gb.

    The disk will then start with size_gb and will resize when the free space
    falls below resize_trigger_gb. The disk will grow by resize_increment_gb
    during each resize, but will not grow beyond max_size_gb.

    Attributes:
        size_gb (int): Disk size in GB.
        max_size_gb (int): Maximum allowed disk size in GB. Disk will not grow
            beyond this limit.
    """
    max_size_gb: int

    # resize_trigger_gb (int): Free space threshold (in GB) to trigger
    #     resizing. Example: If set to 5GB, the disk will resize when 5GB of
    #     free space remains.
    resize_trigger_gb: int = field(default=5, init=False)

    # resize_increment_gb (int): Amount (in GB) to increase the disk size
    #     during each resize.
    resize_increment_gb: int = field(default=10, init=False)
    is_resizable: bool = True

    @property
    def resize_config(self):
        """Return the disk resize configuration as a dictionary.
        If the disk is not resizable, it will return None."""
        if not self.is_resizable:
            return None
        return {
            "free_space_threshold_gb": self.resize_trigger_gb,
            "size_increment_gb": self.resize_increment_gb,
            "max_disk_size_gb": self.max_size_gb,
        }
