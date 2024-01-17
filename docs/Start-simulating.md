# Start simulating

After having installed and prepared your environment, you can start simulating with Inductiva API.

Here, we will show you the first example you can run to test if your setup is working properly. We will use the [Wind Tunnel](Wind-Tunnel) scenario to simulate the flow around a vehicle.

## Example
To run this simulation and render the visualizations install the fluids extra dependencies with

```
!pip install --upgrade inductiva[fluids_extra]
```

Follow the example with:

```python
import inductiva

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
    object_path=vehicle_path,
    num_iterations=50,
    resolution="low")

# Download the simulation output to your local machine.
output = task.get_output()
```

With this code snippet, you will initialise your simulation, submit it to the API, wait for it to finish and download the output files to your local machine.

Then, in some particular cases, you have built-in methods to render the output files. For example, in the case of the Wind Tunnel scenario, you can render the pressure field over the object with:

```
# Get the pressure field over the object
pressure_field = output.get_object_pressure_field()
```

To render the post-processing on a computer with the physical display run:
```python
pressure_field.render()
```

Or, in Colab or a headless server, install the extra dependencies with

```
!apt install libgl1-mesa-glx xvfb
``` 
and run:

```python
pressure_field.render(virtual_display=True, save_path="pressure_field.png")
```

<p align="center">
  <img src="https://github.com/inductiva/inductiva/blob/main/assets/media/openfoam/pressure_field.png?raw=true" alt="Pressure Field of a vehicle." width="400" height="330" >
</p>