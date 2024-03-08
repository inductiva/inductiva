"""Inductiva template engine class to manipulate files and render templates."""

from typing import Optional
import pathlib

from inductiva import types
from inductiva.utils import file_manager
from inductiva.templating import renderers


class TemplateEngine(file_manager.FileManager):
    """Template Engine for managing files and rendering templates.
    
    This class inherits the properties of the file_manager to manage and
    add files and directories. Users can use different rendering engines in
    a common framework to generalize their required files.
    """

    def __init__(self,
                 template_dir: types.Path,
                 root_dir: types.Path = "rendered_dir",
                 renderer: renderers.BaseRenderer = renderers.JinjaRenderer):
        """Initialize the template engine.
        
        Args:
            template_dir (types.Path): Path to the template directory.
            root_dir (types.Path): Path to the root directory.
            renderer (renderers.BaseRenderer): Template Renderer object.
        """

        self.set_root_dir(root_dir)
        self.__root_dir = self.get_root_dir()
        self.__template_dir = pathlib.Path(template_dir)
        self.renderer = renderer

    def render_file(self,
                    source_file: types.Path,
                    target_file: Optional[types.Path] = None,
                    overwrite: bool = False,
                    **render_args):
        """Render a file from the template_dir to the root_dir. 

        Render a source_file from the root template dir to the root dir. The
        target_file is optional, and if not given, it is named as the source
        file without the rendered extension.

        In case, the name of the target_file being rendered already exists, the
        method by default raises a FileExistsError. To overwrite it anyway, set
        the flag `overwrite` to True.

        Args:
            source_file (types.Path): Path to the source file.
            target_file (types.Path): Path to the target file.
            overwrite (bool): Overwrite the target file if it already exists.
            render_args: Arguments to render the template file.
        """

        if not self.renderer.is_template(source_file):
            raise ValueError(f"{source_file} is not a template file.")
        source_file = self.__template_dir / source_file

        if target_file is None:
            target_file = pathlib.Path(source_file).name

        target_file = self.renderer.strip_extension(target_file)
        target_file = self.__root_dir / target_file

        if not overwrite and target_file.exists():
            raise FileExistsError(f"{target_file} already exists.")

        (target_file.parent).mkdir(parents=True, exist_ok=True)
        self.renderer.render_file(target_file=target_file,
                                  source_file=source_file,
                                  **render_args)

        return target_file

    def render_dir(self,
                   source_dir: Optional[types.Path] = None,
                   target_dir: Optional[types.Path] = None,
                   overwrite: bool = False,
                   **render_args):
        """Render a directory to the root_dir and render all template files.
        
        Render the directory `source_dir` in the template directory to the
        `target_dir` inside the root directory from the given `render_args`.
        It keeps the structure from the source directory. 
        
        If `source_dir` is None, the root template_dir is used. If `target_dir`
        is None, the source_dir is rendered into the root directory. All
        template files will have the template suffix stripped from their name
        upon copying.

        In case, the names of the files being rendered to the target_dir already
        exist, method by default raises a FileExistsError. To overwrite them,
        set the flag `overwrite` to True.

        Args:
            source_dir (types.Path): Path to the source directory.
            target_dir (types.Path): Path to the target directory.
            overwrite (bool): Overwrite the target directory if  already exists.
            render_args: Arguments to render the template files.
        """

        if source_dir is None:
            source_dir = self.__template_dir

        if target_dir is None:
            target_dir = "."

        target_dir = self.__root_dir / target_dir

        if not overwrite:
            try:
                self._validate_destination(source_dir, target_dir)
            except FileExistsError as e:
                msg = f"{e}; set `overwrite=True` to overwrite existing files."
                raise FileExistsError(msg) from e

        self.renderer.render_dir(source_dir=source_dir,
                                 target_dir=target_dir,
                                 **render_args)

        return target_dir

    def _validate_destination(self, source_dir: pathlib.Path,
                              target_dir: pathlib.Path):
        """Check if the target directory exists and is empty."""

        source_map = self.renderer.get_template_files(source_dir)

        for subdir, contents in source_map.items():
            target_subdir = target_dir / subdir
            subdir_files = contents["files"] + \
                contents[self.renderer.TEMPLATE_EXTENSION]

            for file in subdir_files:
                target_file = target_subdir / file
                target_name = self.renderer.strip_extension(target_file)
                if target_name.exists():
                    raise FileExistsError(
                        f"File '{target_name}' already exists.")
