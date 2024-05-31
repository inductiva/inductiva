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
class OutputInfo:
    """Information about the contents of a Task output archive."""
    files: List[FileInfo]

    @property
    def n_files(self) -> int:
        return len(self.files)

    @property
    def total_size(self) -> int:
        return sum(file.size for file in self.files)

    @property
    def total_compressed_size(self) -> int:
        return sum(file.compressed_size for file in self.files)

    def __str__(self) -> str:
        size_header = "Size"
        compressed_header = "Compressed"
        name_header = "Name"

        lines = [
            f"Total size: {format_utils.bytes_formatter(self.total_size)}",
            ("Total compressed size: "
             f"{format_utils.bytes_formatter(self.total_compressed_size)}"),
            f"Number of files: {self.n_files}",
            "Contents:",
            f"  {size_header:<12} {compressed_header:<12} {name_header:<12}",
        ]
        for file in self.files:
            lines.append(
                f"  {format_utils.bytes_formatter(file.size):<12}"
                f" {format_utils.bytes_formatter(file.compressed_size):<12}"
                f" {file.name:<12}")

        return "\n".join(lines)
