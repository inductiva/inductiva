# Install

Inductiva package is simple to install, just run on your terminal:

```
pip install --upgrade inductiva
```

This will provide the core functionalities of the API, which allows you to submit jobs, control machines and run simulations. To use the visualization and post-processing tools, you need to install additional optional dependencies specific to different scientific domains: `molecules_extra`, `fluids_extra` or `coastal_extra`. For example, for fluid dynamics:

```
pip install --upgrade "inductiva[fluids_extra]"
```

For any issues with the installation follow the next steps.

## Installation troubleshooting

### Why can't I install the optional packages?
Depending on your shell, you may encounter issues when trying to install optional packages such as `inductiva[molecules_extra]`. This is because certain shells interpret brackets, like those in `[molecules_extra]`, in a special way. To prevent any misinterpretation or errors, enclose the package name and its extras in double quotes. To ensure a successful installation, please use the following command:

```bash
pip install --upgrade "inductiva[molecules_extra]"
```

### Why can't I install Inductiva package?
If installing the package fails, you can retry it on a new Python virtual environment. A [virtual environment](https://docs.python.org/3/library/venv.html) allows you to have a fresh Python environment with isolated dependencies. In your shell, run:

```
python -m venv <venv>
```

In that command, you should replace `<venv>` with the path (*e.g.*, `.venv`) in which you would like to create the environment. Then, to activate the environment (again, correctly replacing `<venv>`), run:

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

After activating the virtual environment, you can install the package as described below:

```
pip install --upgrade inductiva
```