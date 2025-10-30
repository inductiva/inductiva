# templating

Inductiva template manager class to manipulate files and render templates.

### *class* TemplateManager

Bases: `object`

Template manager for rendering files and directories.

Template files contains variables that are replaced by
values provided as keyword arguments to the render methods.
The underlying template engine is the Jinja2 library.

#### *classmethod* render_dir(source_dir: str, target_dir: str, overwrite: bool = False, \*\*render_args)

Renders the entire template directory to the given target_dir,
using the provided rendering args, preserving the directory structure of
the source_dir. Files that are not templates are copied as is.

When target_dir does not exist, it will be created. If any destination
file inside target_dir already exists, a FileExistsError will be
raised when overwrite is False (default). Otherwise, files are
overwritten.

* **Parameters:**
  * **source_dir** (*str*) – Path to the source directory,
    containing template files (and potentially static files to be
    copied without modification).
  * **target_dir** (*str*) – Path to the target directory.
  * **overwrite** (*bool*) – If True, the destination folder will first
    be deleted if it already exists.
    If False (default), a FileExistsError will be raised if
    any file already exists in the target directory.
  * **render_args** – Keyword arguments to render the template files.
* **Raises:**
  **FileExistsError** – If the destination file already exists and
      overwrite is False.

::docsbannersmall
::
