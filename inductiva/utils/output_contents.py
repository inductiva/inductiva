"""Information about the contents of a Task output archive."""
from dataclasses import dataclass
from typing import List

from inductiva.utils import format_utils


@dataclass
class FileInfo:
    """Information about a file in an output archive."""
    name: str
    size: int
    compressed_size: int


@dataclass
class OutputContents:
    """Information about the contents of a Task output archive."""
    size: int
    contents: List[FileInfo]

    @property
    def num_files(self) -> int:
        return len(self.contents)

    def __str__(self) -> str:
        size_header = "Size"
        compressed_header = "Compressed"
        name_header = "Name"

        lines = [
            f"Archive size: {format_utils.bytes_formatter(self.size)}",
            "Contents:",
            f"  {size_header:<12} {compressed_header:<12} {name_header:<12}",
        ]
        for file in self.contents:
            lines.append(
                f"  {format_utils.bytes_formatter(file.size):<12}"
                f" {format_utils.bytes_formatter(file.compressed_size):<12}"
                f" {file.name:<12}")

        return "\n".join(lines)
