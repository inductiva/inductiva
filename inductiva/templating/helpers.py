"""Helper functions for the templating module."""
import os

TEMPLATE_EXTENSION = ".jinja"


def is_template(file: str) -> bool:
    """Check if the given file is of template type."""
    return file.endswith(TEMPLATE_EXTENSION)


def strip_extension(file: str) -> str:
    """Strip the template extension if the given file is a template."""
    if not is_template(file):
        return file
    return os.path.splitext(file)[0]


def get_dir_structure(template_dir: str):
    """
    Generate a dictionary with the structure of the template directory.
    The keys are the subdirectories of the template directory, relative to
    the template directory itself. The values are dictionaries with the keys
    "files" and "templates" containing the lists of the files and templates
    in each subdirectory, respectively.
    Example:
    >>> renderer.get_dir_structure("template_dir")
    {
        ".": {
            "files": ["file0.1"],
            "templates": ["template0.1"]
        },
        "subdir1": {
            "files": ["file1.1", "file1.2"],
            "templates": ["template1.1", "template1.2"]
        }
        ...
    }
    """
    return {
        os.path.relpath(root, start=template_dir): {
            "templates": [f for f in files if is_template(f)],
            "files": [f for f in files if not is_template(f)],
        } for root, _, files in os.walk(template_dir, topdown=True)
    }
