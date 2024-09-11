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
    resize_trigger_gb, resize_increment_gb, and max_size_gb.

    The disk will then start with size_gb and will resize when the free space
    falls below resize_trigger_gb. The disk will grow by resize_increment_gb
    during each resize, but will not grow beyond max_size_gb.

    Attributes:
        size_gb (int): Disk size in GB.
        resize_trigger_gb (int): Free space threshold (in GB) to trigger
            resizing. Example: If set to 5GB, the disk will resize when 5GB of
            free space remains.
        resize_increment_gb (int): Amount (in GB) to increase the disk size
            during each resize.
        max_size_gb (int): Maximum allowed disk size in GB. Disk will not grow
            beyond this limit.
    """
    size_gb: int
    resize_trigger_gb: int = field(default=None)
    resize_increment_gb: int = field(default=None)
    max_size_gb: int = field(default=None)

    def __post_init__(self):
        args = [
            self.resize_trigger_gb, self.resize_increment_gb, self.max_size_gb
        ]
        if any(arg is not None for arg in args) and not all(arg is not None
                                                            for arg in args):
            raise ValueError(
                "If you provide one of resize_trigger_gb, resize_increment_gb, "
                "or max_size_gb, all must be provided.")

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

    @property
    def is_resizable(self):
        """Check if the disk is resizable.
        A disk is resizable if all the resize attributes are provided.
        (resize_trigger_gb, resize_increment_gb, and max_size_gb)"""

        return all(arg is not None for arg in [
            self.resize_trigger_gb, self.resize_increment_gb, self.max_size_gb
        ])
