
[![Python package](https://github.com/inductiva/inductiva/actions/workflows/python-package.yml/badge.svg)](https://github.com/inductiva/inductiva/actions/workflows/python-package.yml)

![linkedin_header](https://user-images.githubusercontent.com/104431973/231184851-0ce34289-593e-4832-aaa2-9aae652113f5.jpg)

# Inductiva API Python client

Inductiva is a Python package for executing large-scale simulations of physical systems directly in the cloud.

Inductiva API offers distinct advantages:

- 🔄 It consolidates various simulation domains, including fluid and molecular dynamics, plasmas, and structural mechanics, under a single unified entry point.
- 📦 Eliminates the need for installing and managing complex simulation software and corresponding dependencies.
- 🚀 Allows running hundreds or even thousands of simulations concurrently, with no coding.
- 💽 Automatically optimizes hardware configurations for each type of simulation (e.g., CPU vs. GPU, appropriate number of CPU cores, RAM, etc.).
- 🐍 With Inductiva API, you are not limited to a pre-defined GUI or intricate configuration languages and scripts. Instead, you write small python programs that seamlessly integrate with your existing codebase and ML framework.


## Simulators

**Inductiva API** has available several open-source simulators ready to use. Users 
familiar with the simulators can easily start running simulations with their 
previously prepared simulation configuration files. In this way, they can take 
advantage of performant hardware to speed up their simulation and exploration.

Check the following pages for additional details on how to use them with
**Inductiva API**:
- [SPlisHSPlasH](https://github.com/inductiva/inductiva/wiki/SPlisHSPlasH)
- [DualSPHysics](https://github.com/inductiva/inductiva/wiki/DualSPHysics)
- [OpenFOAM](https://github.com/inductiva/inductiva/wiki/OpenFOAM)
- [SWASH](https://github.com/inductiva/inductiva/wiki/SWASH)
- [XBeach](https://github.com/inductiva/inductiva/wiki/XBeach)
- [Reef3D](https://github.com/inductiva/inductiva/wiki/Reef3D)
- [GROMACS](https://github.com/inductiva/inductiva/wiki/GROMACS)
- [FDS](https://github.com/inductiva/inductiva/wiki/FDS)

To learn how to use these simulators with Inductiva API, check the example below and the [Simulators section](https://github.com/inductiva/inductiva/wiki/Simulators).

If you would like other simulators to be added, contact us at [simulations@inductiva.ai](mailto:simulations@inductiva.ai).

### Example

Example of how to use the simulators:

```python
import inductiva

# Download the configuration files
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-input-example.zip", unzip=True
)

# Initialize the Simulator
simulator = inductiva.simulators.REEF3D()

# Run the simulation
task = simulator.run(input_dir=input_dir)
```

The user must specify the input directory containing the files to run the simulation. In the above example, a directory with the configuration of a simulation is downloaded, and passed as argument to the simulator call.


## Installation

Inductiva package is simple to install, just run on your terminal:

```
pip install --upgrade inductiva
```

This will provide the core functionalities of the API, which allows you to submit jobs, control machines and run simulations. 

If you have issues with the installation, check the [Installation troubleshooting](#installation-troubleshooting) for more information.

## API access tokens

Please [request API token](https://docs.google.com/forms/d/e/1FAIpQLSflytIIwzaBE_ZzoRloVm3uTo1OQCH6Cqhw3bhFVnC61s7Wmw/viewform) and set your
API key as an environment variable in your terminal as follows:

```bash
export INDUCTIVA_API_KEY="YOUR_API_KEY"
```

And you are good to go! You can start exploring Inductiva API with the examples below.

## More info:

- [Managing submitted tasks](https://github.com/inductiva/inductiva/tree/main/inductiva/tasks#tasks)
- [Managing computation resources](https://github.com/inductiva/inductiva/tree/main/inductiva/resources#manage-computational-resources)

## Installation troubleshooting

### Why can't I install Inductiva package?
If installing the package failed, you can retry it on a new Python virtual environment. A [virtual environment](https://docs.python.org/3/library/venv.html) allows you to have a fresh Python environment with isolated dependencies. In your shell, run:

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
