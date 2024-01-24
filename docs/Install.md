# Installing Inductiva

Inductiva package is simple to install, just run on your terminal:

```
pip install --upgrade inductiva
```

This will provide the core functionalities of the API, which allows you to submit 
jobs, control machines, and run simulations.

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
