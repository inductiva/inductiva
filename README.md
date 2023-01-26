[![Python package](https://github.com/{{GITHUB_REPOSITORY}}/actions/workflows/python-package.yml/badge.svg)](https://github.com/{{GITHUB_REPOSITORY}}/actions/workflows/python-package.yml)

[![Documentation Status](https://readthedocs.com/projects/inductiva-research-labs-{{GITHUB_REPOSITORY}}/badge/?version=latest&token=a7a7344c5bc9e511ca715449f94612e4bb7c2dff4c90931778b9c85dbbf5480c)](https://inductiva-research-labs-{{GITHUB_REPOSITORY}}.readthedocs-hosted.com/en/latest/?badge=latest)

# Start a project

Use this template to generate a repository for your new project. By doing so you
will quickly have access to a couple of tools we've already set up for you.

## Python

We love Python! To make our lives easier and write good code we follow a couple
of common best practices:

* For local development we recommend using Python's virtual environments. This
  way you can keep a separate environment for each project without polluting
  your operating system with additional Python packages. We also recommend using
  `pip` for managing those packages.

  To create a new virtual environment (only once):

  ```bash
  python3 -m venv .env
  ```

  To activate your virtual environment (every time you want to work on the
  current project):

  ```bash
  source .env/bin/activate 
  ```

  It's always a good practice to upgrade `pip` when you create the virtual environment, before installing other packages. In the past we've noticed some
  bugs with older versions of `pip`, causing the installation of incompatible
  versions of several packages.

  ```bash
  pip install --upgrade pip
  ```

  To install the project's required packages (whenever dependencies must be
  updated):

  ```bash
  pip install -r requirements.txt
  ```

* We follow
  [Google's style guide](https://google.github.io/styleguide/pyguide.html).
  Specifically, we format code with `yapf`, and check for bugs and style
  problems with `pylint`.

  Usually, there's no need to run these tools manually. We added continuous
  integration (CI) through GitHub, so every time your pull request has errors
  you'll be alerted via email. Over time you'll get used to these rules and
  coding consistently will become second nature.

  However, if you want to run `yapf` and `pylint` locally, simply install them
  via `pip`:

  ```bash
  pip install yapf pylint
  ```

* We write tests using `pytest`. See [`example_test.py`](example_test.py) for an
  example and delete if example after you copied the template to your new
  project.

  Tests are also run automatically in GitHub via continuous integration. If you
  want to run them locally, install `pytest`:

  ```bash
  pip install pytest
  ```

  and run:

  ```bash
  pytest
  ```

  in your base directory.

  `pytest` is smart enough to discover your tests automatically when you're
  following
  [the right conventions](https://docs.pytest.org/en/stable/goodpractices.html#conventions-for-python-test-discovery). At Inductiva we'd prefer that any
  test file (e.g. `example_test.py`) stays in the same directory next to the
  component it tests (e.g. `example.py`).

## Documentation

We use [Sphinx](https://www.sphinx-doc.org/) to generate documentation to our
projects. The documentation files can be written using either Markdown or
reStructuredText syntax and are stored in the `docs` directory. See
[`docs/index.md`](docs/index.md) for a collection of simple examples to get you
started including math formulas and code snippets in your documentation.

If you want to generate the documentation locally, simply install the
dependencies via `pip`:

```bash
pip install -r docs/requirements.txt
```

and run

```bash
sphinx-build -b html docs docs/_build
```

in your base directory to build the documentation in `docs/_build`.

To serve the documentation locally and automatically build it with new changes,
install the `sphinx-autobuild` package,

```bash
pip install sphinx-autobuild
```

and run in your base directory

```bash
sphinx-autobuild -b html docs docs/_build
```

This will serve the documentation at <http://127.0.0.1:8000> and automatically
update it when you save changes to the documentation files.

The documentation may also be made available on
[Read The Docs](https://readthedocs.org/) for configured repositories. When set
up, this integration will build and deploy the documentation on each commit to
the `main` branch.

See the documentation for the start-a-project repo at:

https://inductiva-research-labs-start-a-project.readthedocs-hosted.com/en/latest/

Apart from documentation files, we also prepare README files using Markdown
syntax to present short introductions of the projects. This README is one such
example.


## Conda Environments

To automate all the setup commands above, we have created a Conda
environment config file that you can use for your convenience.

Make sure the first line of the conda.yml file contains a name related to 
your project, such as:

```
name: startaproject_env
```

Note: If you are starting a new project from this template repo, it might work 
automatically for you, because we use the placeholder `name: {{GITHUB_REPOSITORY}}_env`
that gets auto-filled with the repo name.

Now, if you have conda installed, just run:

```bash
conda env create -f conda.yml
```

This will install python 3.9, pylint, yapf, pytest, sphinx, and all the 
package requirements, and you can now activate that environment with:

```bash
conda activate startaproject_env
```

And you are good to go!