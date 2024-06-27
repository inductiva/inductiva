# Troubleshooting Guide
In this troubleshooting section, we aim to guide you through resolving issues
related to incorrect or incomplete installation of the Inductiva package
or its dependencies. This guide assumes you’ve gone through all steps we've shared
in the [user console](https://genesis.inductiva.ai), and have faced
some issues.

If you find bugs, need help, or want to talk to the developers, reach out to us on
support@inductiva.ai

If you find any security issues, please report to security@inductiva.ai

## Installation Failures
If you encounter issues installing the Inductiva API package, there are several
steps you can take:

### Ensure you have a working pip
```
pip install --upgrade pip
```
### Use virtualenv or venv to isolate dependencies

If installing the package fails, you can retry it on a new Python virtual environment.
A [virtual environment](https://docs.python.org/3/library/venv.html) allows you to
have a fresh Python environment with isolated dependencies.

In your shell, run:

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
