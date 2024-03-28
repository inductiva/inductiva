# Templating Manager

## Introduction

The Inductiva API is all about enabling you to simulate at scale. As we have shown,
with a few lines of Python code, you can send your simulations to MPI Clusters
assembled from last-generation cloud hardware, letting you run much larger simulations
than you would be able to use your local resources. Or you can spin up a large Machine
Group, with dozens or hundreds of machines, and send a large number of simulations
to be run on those machines in parallel. Such massive parallelism is precisely
what you need when you want to find the optimal solution for a problem, and the
only way to test each candidate's solution is by simulating it. To do so, one needs
to configure multiple simulations, each with a slightly different set of parameters.
Generally speaking, the goal is to explore the largest possible extent of the
design space.

In this context, the Inductiva API provides a powerful tool for exploring the design
space of a problem: the *templating manager*. The templating manager allows you to
quickly generate a large number of simulation configurations by starting from a
base case and replacing some of its fixed values with variables that you can programmatically
change.

## The `TemplateManager` class

The `TemplateManager` class is a utility class that allows you to manage templating
files and specify how to render each template into a concrete configuration file.
It abstracts the rendering process, allowing the user to focus solely on defining
the source of the template files, the destination directory of the rendered files,
and the values of the variables to be used in the rendering process.

In the following sections, we'll provide a basic introduction to rendering concepts
and explain how to use the `TemplateManager` class to render different variations
of template files to a destination folder. In the end, we will discuss some safety
features that the `TemplateManager` class provides to prevent accidental overwriting
of files and ensure the uniqueness of the destination directory.

### Rendering basics

To start, a template file is a file that contains labels that will be replaced
with specific values during a rendering process. These labels -- variables --
are enclosed in double curly brackets, optionally specifying default values they
might take when not explicitly set during the rendering process.

In the following example, the content of a template file defines a configuration
parameter named `config_parameter`:

```jinja
config_parameter = {{ parameter_value }}
```

When this file is rendered with `10` for the variable `parameter_value`,
the resulting file will read:

```txt
config_parameter = 10
```

If the template file were to be defined with a default value for `parameter_value`:

```jinja
config_parameter = {{ parameter_value | default(20) }}
```

and **rendered without providing a value for it**, the resulting file would read:

```txt
config_parameter = 20
```

### Rendering with the `TemplateManager` class

The `TemplateManager` class is initialized with the path to the directory containing
the template files (the *template_dir*) and the name of the directory where the
rendered files will be stored (the *root_dir*). The latter is optional and defaults
to a directory named `rendered_dir` inside the current working directory.
If the output directory does not exist, it will be created.

```python
import inductiva

# New instance of the TemplateManager class specifying the name 
# of the template directory. All rendered files will be stored
# inside ./rendered_dir (the default root directory of the manager).
template_manager = inductiva.TemplateManager(template_dir="my_templates")

# new instance of the TemplateManager class specifying both the name 
# of the template directory and the name of the destination directory
# of the rendered files. In this case, all rendered files will be stored
# inside ./my_rendered_files.
template_manager = inductiva.TemplateManager(template_dir="my_templates",
                                             root_dir="my_rendered_files")
```

The `TemplateManager` class exposes two rendering methods: `render_dir` and
`render_file`. The former renders all files inside the template directory,
while the latter only renders a single file. Both methods take a set of keyword
arguments that specify the values that the variables in the template files will
take. The manager identifies template files by looking for the `.jinja` extension
in the file name. Upon rendering, the extension is removed from the file name inside
the destination directory.

Let's now look into details on how to use these two methods:

### Rendering a directory

For the sake of illustration, let's consider a directory, containing three template
files and a regular (non-template) file, with the following structure:

```console
$ tree my_templates
my_templates
├── config.yaml.jinja
├── readme.txt
└── objects
    ├── obj1.def.jinja
    └── obj2.def.jinja
```

The contents of each file are as follows:

```console
$ cat my_templates/config.yaml.jinja
density: {{ density }}
viscosity: {{ viscosity }}
position: {{ position | default([0, 0, 0]) }}

$ cat my_templates/readme.txt
I am a non-template file. I will be copied as is.

$ cat my_templates/objects/obj1.def.jinja
type: sphere
radius: {{ sphere_radius }}

$ cat my_templates/objects/obj2.def.jinja
type: cube
side: {{ cube_side | default(1) }}
```

To render all files in the `my_templates` directory, you can use the `render_dir`:

```python
import inductiva

# instantiate the TemplateManager, specifying the template directory
# and the name of the destination directory
template_manager = inductiva.TemplateManager(template_dir="my_templates",
                                             root_dir="my_rendered_files")

# render all files in the template directory specifying the values of
# the variables in the template files. Note that we are deliberately not
# providing values for the `position` or `cube_side` variables in the
# config.yaml.jinja and obj2.def.jinja files, respectively. This
# will enfore the default values to be used.
template_manager.render_dir(density=1000, viscosity=1e-6, radius=0.5)
```

After running the code above, the `my_rendered_files` directory will contain the
following files:

```console
$ tree my_rendered_files
my_rendered_files
├── config.yaml
├── readme.txt
└── objects
    ├── obj1.def
    └── obj2.def

$ cat my_rendered_files/config.yaml
density: 1000
viscosity: 1e-6
position: [0, 0, 0]

$ cat my_templates/readme.txt
I am a non-template file. I will be copied as is.

$ cat my_rendered_files/objects/obj1.def
type: sphere
radius: 0.5

$ cat my_rendered_files/objects/obj2.def
type: cube
side: 1
```

One can optionally request the rendering of a subdirectory inside the template
directory by providing the name of the source subdirectory as an argument to the
`render_dir` method. At the same time, one can specify the name of the destination
folder inside the root output directory.
This mechanism is useful when the template directory contains template files
for different simulation scenarios/stages but you want to keep the flexibility
of selecting which scenario/stage to render. For example, in the following snippet,
we only render the `objects` subdirectory to `my_rendered_files/only_objects`:

```python
template_manager = ...
template_manager.render_dir('objects', # subdirectory to render
                            'only_objects', # destination subdirectory
                            radius=0.5, side=2)
```

```console
$ tree my_rendered_files
my_rendered_files
└── only_objects
    ├── obj1.def
    └── obj2.def
```

### Rendering a single file

Sometimes, you may just want to render a single file from the template directory.
In this case, you can use the `render_file` method. This method takes the name of
the file to render as an argument along with the values for the variables in the
template file.

Using the same template files as before, let's render the `obj1.def.jinja` file.
By default, the rendered file will be saved to the root output directory:

```python
template_manager = ...
template_manager.render_file('objects/obj1.def.jinja',
                            radius=10)
```

```console
$ tree my_rendered_files
└── obj1.def.yaml
```

Similarly to the `render_dir` method, you can specify the destination name
of the rendered file by providing the `target_file` argument:

```python
template_manager = ...
template_manager.render_file('objects/obj1.def.jinja',
                             target_file='my_objects/sphere.def',
                             radius=10)
```

```console
$ tree my_rendered_files
my_rendered_files
└── my_objects
    └── sphere.def
```

### Adding external resources

In addition to rendering template files, the `TemplateManager` class also provides
mechanisms to add external files to the destination directory. This method is useful
when you need to copy files that are not part of the template directory but are
required for the simulation. The `copy_dir` and `copy_file` methods allow you to
copy entire directories or individual files to the destination directory, respectively.

In the following examples, we add external resources to the destination directory,
first by adding an entire directory and then by adding a single file:

```python
template_manager = ...
template_manager.copy_dir('/path/to/external_resources')
template_manager.copy_file('./external_file.txt')
```

```console
$ tree my_rendered_files
my_rendered_files
├── external_file.txt
└── external_resources
    ├── ...
```

Similarly to the `render_dir` and `render_file` methods, you can specify the
destination name of the copied directory or file:

```python
template_manager = ...
template_manager.copy_dir('/path/to/external_resources',
                          'resources/copied_resources')
template_manager.copy_file('./external_file.txt',
                           'copied_file.txt')
```

```console
$ tree my_rendered_files
my_rendered_files
├── copied_file.txt
└── resources
    └── copied_resources
        ├── ...
```

### Overwrite safety

By default, the `TemplateManager` *will not overwrite* any existing files in the
destination directory. Calls to the `render_*` or `copy_*` methods will fail if
any of the destination files already exist. Rendering and coping actions are
transactional, meaning that the entire action will fail if any destination file
exists. For example, consecutive calls to the same method will fail in the second
call. This behavior is intended to prevent accidental overwriting of
files that may have been generated in a previous run.

To enforce the overwriting of existing files, you can set the `overwrite` argument
to `True` when calling the `render_*` or `copy_*` methods.

```python
template_manager = ...

# ✔ this call will succeed because the destination directory is empty
template_manager.render_dir(density=1000, viscosity=1e-6,
                            position=[0, 0, 0], radius=0.5)

# ✖ the second call will fail with an FileExistsError because at least
# one of the rendered files would overwrite the equivalent file
# generated in the first call
template_manager.render_dir(density=1000, viscosity=1e-6,
                            position=[0, 0, 0], radius=0.5)
                        
# ✔ to ensure the call succeeds, set the `overwrite` argument to `True`
template_manager.render_dir(overwrite=True,
                            density=1000, viscosity=1e-6,
                            position=[0, 0, 0], radius=0.5)
```

### Uniqueness of the destination directory

The `TemplateManager` class ensures that the destination directory is unique for
each instance. This means that if you instantiate two `TemplateManager` objects
with the same destination directory, the second object will point to a slightly
different destination directory.

```python
>>> manager1 = inductiva.TemplateManager(..., root_dir="my_rendered_files")
>>> manager2 = inductiva.TemplateManager(..., root_dir="my_rendered_files")
>>> print(manager1.get_root_dir())
my_rendered_files
>>> print(manager2.get_root_dir())
my_rendered_files__2
```

This way, the `TemplateManager` guarantees that the destination directory is unique
and that no files are accidentally overwritten. This behavior can be changed by
setting the `INDUCTIVA_DISABLE_FILEMANAGER_AUTOSUFFIX` environment variable to
`True` before instantiating the `TemplateManager` object. In this case, if the
destination directory already exists, a `FileExistsError` exception will be thrown
when instantiating the `TemplateManager` object.

```python
>>> os.setenv('INDUCTIVA_DISABLE_FILEMANAGER_AUTOSUFFIX', 'True')
>>> manager1 = inductiva.TemplateManager(..., root_dir="my_rendered_files")
>>> manager2 = inductiva.TemplateManager(..., root_dir="my_rendered_files")
---------------------------------------------------------------------------
FileExistsError                           Traceback (most recent call last)
line 1
----> 1 raise FileExistsError(f"Directory {root_dir} already exists.")

FileExistsError: Directory my_rendered_files already exists.
```

When using the templating manager inside a loop, it is important to ensure that
the destination directory is unique for each iteration. This can be achieved by
setting the `root_dir` argument to a unique value for each iteration or by relying
on the above mechanism to ensure uniqueness across different iterations.

```python
template_manager = inductiva.TemplateManager(..., root_dir=)

# explicitly set an unique root directory for each iteration
for iteration in range(...):
    template_manager.set_root_dir(f"my_rendered_files_iter{iteration}")
    print(template_manager.get_root_dir())

# would print:
# my_rendered_files_iter0
# my_rendered_files_iter1
# my_rendered_files_iter2
# ....

# or let the manager define an unique root directory each iteration
# (assuming the INDUCTIVA_DISABLE_FILEMANAGER_AUTOSUFFIX environment
# variable is either undefined or set to false)
for iteration in range(...):
    template_manager.set_root_dir("my_rendered_files")
    print(template_manager.get_root_dir())

# would print:
# my_rendered_files
# my_rendered_files__2
# my_rendered_files__3
# ....

```
