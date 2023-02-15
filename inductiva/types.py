"""Simple classes used for type definitions."""
import os

class DirPath():
    """Util class to represent paths to directories.

    When calling API methods, the type annotation for this class is used to
    detect that the input is a directory, so that the `pack_input` function
    can handle it correctly.

    Properties:
        path: String that represents the path to the directory.
    """
    def __init__(self, path: str):
        if not os.path.isdir(path):
            raise ValueError("The provided path is not a directory.")
        self.path = os.path(path)
