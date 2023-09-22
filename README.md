
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


## Installation

Inductiva package is simple to install, just run on your terminal:

```
pip install --upgrade inductiva
```

This will provide the core functionalities of the API, which allows you to submit jobs, control machines and run simulations. To use the visualization and post-processing tools, you need to install additional optional dependencies specific to different scientific domains: `molecules_extra`, `fluids_extra` or `coastal_extra`. For example, for fluid dynamics:

```
pip install --upgrade "inductiva[fluids_extra]"
```

If you have issues with the installation, check the [Installation troubleshooting](#installation-troubleshooting) for more information.

## API access tokens

Please [request API token](https://docs.google.com/forms/d/e/1FAIpQLSflytIIwzaBE_ZzoRloVm3uTo1OQCH6Cqhw3bhFVnC61s7Wmw/viewform) and add the following line to your code:

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"
```

And you are good to go! You can start exploring Inductiva API with the examples below.

## Pre-built Simulation Scenarios

**Inductiva API** contains pre-built simulation scenarios that define physical systems of interest ready to simulate. Users can choose some parameters and configure the system according to their needs, run the simulation using the most adequate resources and visualize the results.

### WindTunnel Example

To run this simulation and render the visualizations install the fluids extra dependencies with
```python
!pip install --upgrade inductiva[fluids_extra]
```

Follow the example with:

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Url to a test object in Inductiva Github repository
vehicle_url = "https://raw.githubusercontent.com/inductiva/inductiva/main" \
              "/resources/vehicle.obj"
vehicle_path = inductiva.utils.files.download_from_url(vehicle_url)

# Initialize the scenario
scenario = inductiva.fluids.WindTunnel(
    flow_velocity=[30, 0, 0],
    domain={"x": [-5, 15], "y": [-5, 5], "z": [0, 8]})

# Run a simulation
task = scenario.simulate(
    object_path=vehicle_path, num_iterations=50, resolution="low")

# Download the simulation output to your local machine.
output = task.get_output()

# Get the pressure field over the object
pressure_field = output.get_object_pressure_field()
```

To render the post-processing on a computer with physical display run:
```python
pressure_field.render()
```

Or, in Colab or headless server, install the extra dependencies with

```
!apt install libgl1-mesa-glx xvfb
``` 
and run:

```python
pressure_field.render(virtual_display=True, save_path="pressure_field.png")
```

<p align="center">
  <img src="https://github.com/inductiva/inductiva/blob/main/resources/media/openfoam/pressure_field.png?raw=true" alt="Pressure Field of a vehicle." width="400" height="330" >
</p>

### Available Simulation Scenarios

These are the currently available scenarios:


<div align="center">
  
| [**Protein Solvation**](https://github.com/inductiva/inductiva/tree/main/inductiva/molecules/protein_solvation#protein-solvation-scenario) | [**Coastal Area**](https://github.com/inductiva/inductiva/tree/main/inductiva/coastal#coastal-area-scenario)  |
| :-------------------------: | :-------------------------------: |
|<div align="center"><img src="https://github.com/inductiva/inductiva/blob/main/resources/media/md/protein_solvation_big_molecule.gif?raw=true" alt="Protein Solvation simulation" width="200" height="200" /></div> |  <div align="center"><img src="https://github.com/inductiva/inductiva/blob/main/resources/media/random_coastal_area.gif?raw=true" alt="Coastal area simulation" width="250" height="150" /></div>  |
| [**Fluid Tank**](https://github.com/inductiva/inductiva/tree/main/inductiva/fluids/fluid_tank#fluid-tank-scenario) | [**Wind Tunnel**](https://github.com/inductiva/inductiva/tree/main/inductiva/fluids/wind_tunnel#wind-tunnel-scenario) |
| <div align="center"><img src="https://github.com/inductiva/inductiva/blob/main/resources/media/fluid_tank.gif?raw=true" alt="Fluid Tank simulation" width="250" height="150" /></div> | <div align="center"><img src="https://github.com/inductiva/inductiva/blob/main/resources/media/openfoam/default_pressure_field.png?raw=true" alt="Wind tunnel simulation" width="300" height="200" /></div> |

</div>


## Simulators

**Inductiva API** has available several open-source simulators ready to use. Users familiar with the simulators can easily start running simulations with their previously prepared simulation configuration files. In this way, they can take advantage of performant hardware to speed up their simulation and exploration.

The simulators we provide are all open-source and have their own dedicated documentation:
- [SPlisHSPlasH](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH)
- [DualSPHysics](https://github.com/DualSPHysics/DualSPHysics)
- [OpenFOAM](https://www.openfoam.com/)
- [SWASH](https://swash.sourceforge.io/)
- [XBeach](https://oss.deltares.nl/web/xbeach/)
- [GROMACS](https://www.gromacs.org/)

To learn how to use these simulators with Inductiva API, check the example below and the [Simulators section](https://github.com/inductiva/inductiva/tree/main/inductiva/simulators#simulators).

If you would like other simulators to be added, contact us at [simulations@inductiva.ai](mailto:simulations@inductiva.ai).

### Example

Example of how to use the simulators:

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Download the configuration files
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "dualsph-flow-cylinder.zip"
)

# Initialize the Simulator
simulator = inductiva.simulators.DualSPHysics()

# Run the simulation
task = simulator.run(input_dir=input_dir)
```

The user must specify the input directory containing the files to run the simulation. In the above example, a directory with the configuration of a simulation is downloaded, and passed as argument to the simulator call.


## Running multiple simulations in parallel

Up until now, all simulations were run synchronously which gives feedback while the simulation is running until it finishes. However, this is not always the best option, for example, when launching two simulations it implies that they run one after the other. This
is not optimal for when users need to launch hundreds of simulations.

To solve this problem **Inductiva API** allows users to run simulations asynchronously. This means that the simulation is launched and the user can continue with other tasks - like launching more simulations. 

Let's look at an example using the wind tunnel scenario:

```python
import inductiva

vehicle_url = "https://raw.githubusercontent.com/inductiva/inductiva/main" \
              "/resources/vehicle.obj"
vehicle_path = inductiva.utils.files.download_from_url(vehicle_url)

task_list = []
velocity_list=[1, 10, 20, 30, 40]
for velocity in velocity_list:
  scenario = inductiva.fluids.WindTunnel(flow_velocity=[velocity, 0, 0])
  task = scenario.simulate(object_path=vehicle_path, run_async=True)
  task_list.append(task)
```

All simulations are launched in one go, allowing users to continue working on other things. To monitor the progress of individual simulations, users can use `task.get_status()`, or they can view a list of the most recent tasks launched by using `inductiva.tasks.list(last_n=5)`.

Finally, to retrieve the results the user can use `task.get_output()`, which waits for the simulation to finish before downloading the results. During this waiting period, it temporarily blocks the execution of other code. Check the [Tasks section](https://github.com/inductiva/inductiva/tree/main/inductiva/tasks#tasks) for more information on how to do this.

## More info:

- [Managing submitted tasks](https://github.com/inductiva/inductiva/tree/main/inductiva/tasks#tasks)
- [Managing computation resources](https://github.com/inductiva/inductiva/tree/main/inductiva/resources#manage-computational-resources)

## Installation troubleshooting
### Why can't I install the optional packages?
Depending on your shell, you may encounter issues when trying to install optional packages such as `inductiva[molecules_extra]`. This is because certain shells interpret brackets, like those in `[molecules_extra]`, in a special way. To prevent any misinterpretation or errors, enclose the package name and its extras in double quotes. To ensure a successful installation, please use the following command:

```bash
pip install --upgrade "inductiva[molecules_extra]"
```

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
