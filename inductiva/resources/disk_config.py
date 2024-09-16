"""
This module contains the DiskConfig class which is used to store the
configuration for the disk.
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class DiskConfig:
    """Dataclass to represent the disk configuration.
    
    A disk has a size in GB and can be resized based on the free space left.
    To enable resizing, the following attributes must be provided:
    free_space_threshold_gb, size_increment_gb, and max_size_gb.

    The disk will then start with size_gb and will resize when the free space
    falls below free_space_threshold_gb. The disk will grow by
    size_increment_gb during each resize, but will not grow beyond
    max_size_gb.

    Attributes:
        size_gb (int): Disk size in GB.
        max_size_gb (int): Maximum allowed disk size in GB. Disk will not grow
            beyond this limit.
    """
    max_size_gb: int
    is_resizable: bool = True

    # free_space_threshold_gb (int): Free space threshold (in GB) to trigger
    #     resizing. Example: If set to 5GB, the disk will resize when 5GB of
    #     free space remains.
    free_space_threshold_gb: int = field(default=5, init=False)

    # size_increment_gb (int): Amount (in GB) to increase the disk size
    #     during each resize.
    size_increment_gb: int = field(default=10, init=False)

    def __post_init__(self):
        if not isinstance(self.max_size_gb, int) or self.max_size_gb < 1:
            raise ValueError("`max_size_gb` must be a positive integer.")
        if not isinstance(self.is_resizable, bool):
            raise ValueError("`is_resizable` must be a boolean.")
