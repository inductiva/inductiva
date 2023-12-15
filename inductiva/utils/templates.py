"""Utils related to template files."""

import os
import pathlib
import glob
import shutil

import jinja2

from inductiva.utils.files import find_path_to_package

TEMPLATES_PATH = find_path_to_package("templates")


def render_from_string(template_string, **render_args):
    """Render commands with the given arguments."""

    template_env = jinja2.Environment(loader=jinja2.BaseLoader())
    template = template_env.from_string(template_string)
    string = template.render(**render_args)

    return string


def render_file(source_file, target_file, remove_template=False, **render_args):

    source_path = pathlib.Path(source_file)

    source_dir = source_path.parent
    source_file = source_path.name

    environment = jinja2.Environment(loader=jinja2.FileSystemLoader(source_dir))
    template = environment.get_template(source_file)
    stream = template.stream(**render_args)
    stream.dump(target_file)

    if remove_template:
        os.remove(source_path)


def render_directory(source_dir,
                     target_dir,
                     remove_template_dir=False,
                     **render_args):

    shutil.copytree(source_dir, target_dir, dirs_exist_ok=True, symlinks=True)

    jinja_files = get_jinja_files(target_dir)

    for jinja_file in jinja_files:
        jinja_path = os.path.join(target_dir, jinja_file)
        target_path = jinja_path.split(".jinja")[0]

        render_file(source_file=jinja_path,
                    target_file=target_path,
                    remove_template=True,
                    **render_args)

    if remove_template_dir:
        shutil.rmtree(source_dir)


def get_jinja_files(src_dir):
    """Get all jinja files in a directory."""

    jinja_paths = glob.glob(os.path.join(src_dir, "**"
                                         "*.jinja"),
                            recursive=True)

    jinja_files = [
        os.path.relpath(file_path, start=src_dir) for file_path in jinja_paths
    ]

    return jinja_files
