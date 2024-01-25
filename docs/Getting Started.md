# Getting Started
some text here to intoduce this section

## Prerequisites
You will need [Python](https://www.python.org/) 3.8 or higher installed. If you do not have [Python](https://www.python.org/) installed already, here are the instructions to [download Python](https://www.python.org/downloads/) and install it.

> **Tip:** New versions of Python are released annually in October, and it can take a few months for the scientific Python ecosystem to catch up. If you have trouble installing plasmapy on the most recent Python version between October and March, then try installing it on the second most recent version.



## Installation

To install the latest Inductiva package release on PyPI with pip, just open your 
terminal and run:

```
pip install --upgrade inductiva
```

This will provide the core functionalities of the API, which allows you to submit 
jobs, control machines, and run simulations.

**You're all set!** Visit the [example gallery]() to get inspired on the many ways you can use Inductiva. It 
takes about 10 minutes to browse.

For any issues with the installation, see [Troubleshooting](#troubleshooting).

## Troubleshooting

If installing the package fails, you can retry it on a new Python virtual environment. 
A [virtual environment](https://docs.python.org/3/library/venv.html) allows you to 
have a fresh Python environment with isolated dependencies. In your shell, run:

```
python -m venv <venv>
```

In that command, you should replace `<venv>` with the path (*e.g.*, `.venv`) in 
which you would like to create the environment. Then, to activate the environment 
(again, correctly replacing `<venv>`), run:

For `bash`/`zsh`:

```
source <venv>/bin/activate
```

For `cmd.exe` (Windows):

```
<venv>\Scripts\activate.bat
```

For `PowerShell` (Windows):
```
<venv>\Scripts\Activate.ps1
```

After activating the virtual environment, you can install the package as described 
below:

```
pip install --upgrade inductiva
```

## Frequently Asked Questions
