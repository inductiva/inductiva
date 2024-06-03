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
class TaskOutputInfo:
    """Information about the contents of a Task output archive."""
    files: List[FileInfo]
    total_size_bytes: int

    @property
    def n_files(self) -> int:
        return len(self.files)

    @property
    def total_compressed_size_bytes(self) -> int:
        return sum(file.compressed_size for file in self.files)

    def __str__(self) -> str:
        size_header = "Size"
        compressed_header = "Compressed"
        name_header = "Name"

        total_size_formatted = format_utils.bytes_formatter(
            self.total_size_bytes)
        total_compressed_size_formatted = format_utils.bytes_formatter(
            self.total_compressed_size_bytes)

        lines = [
            f"Total size: {total_size_formatted}",
            f"Total compressed size: {total_compressed_size_formatted}",
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
