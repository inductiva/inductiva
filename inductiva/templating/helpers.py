"""Helper functions for the templating module."""
import os
import pathlib

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


def _check_prerender_dir(source_dir_struct, target_dir: str):
    """Check if the destination filenames exist.

    Check if the destination filenames exist and raise a FileExistsError if
    any of them already exists. For template files, the extension is
    stripped before checking.

    Args:
        source_dir_struct (dict): The structure of the source directory,
            as returned by `get_dir_structure`.
        target_dir (pathlib.Path): The destination directory.
    """
    target_dir = pathlib.Path(target_dir)
    for subdir, contents in source_dir_struct.items():
        dest_subdir = target_dir / subdir
        for file in contents["files"] + contents["templates"]:
            raw_target_name = str(dest_subdir / file)
            target_name = strip_extension(raw_target_name)
            if os.path.exists(target_name):
                raise FileExistsError(f"File {target_name} already exists.")
