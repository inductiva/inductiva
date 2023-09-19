# Wind Tunnel Scenario

This scenario models the aerodynamics of an object inside a virtual wind tunnel for a given airflow velocity. Air is injected through one side of the wind tunnel and exits through an outlet on the opposite side. The air flow within the tunnel is modified according to the structure of the object. The system is modelled with the steady-state equations for incompressible flow and the $k-\epsilon$ turbulence models. The simulation is currently performed with the OpenFOAM simulator.

To initialize the scenario, the user can define the following parameters:
- `domain`: Dimensions of the wind tunnel in meters, e.g. `{"x": [-5, 15], "y": [-5, 5], "z": [0, 8]}`.
- `flow_velocity`: Airflow velocity in m/s, e.g. `[30, 0, 0]`. Notice that the airflow is injected from the negative x-direction.

Now, the user is ready to simulate the steady state. Here, the user chooses the object to be inserted inside the wind tunnel. This object is defined with a geometry file in [STL](https://en.wikipedia.org/wiki/STL_(file_format)) or [OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) format. Notice that:
- the object is not re-scaled automatically, so the user must either change the dimensions of the domain or scale the object himself.
- the required meshing step is automatically done on the backend before the simulation starts with [snappyHexMesh tool](https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.4-mesh-generation-with-the-snappyhexmesh-utility) of OpenFOAM.

To test this scenario we make available one example, but users can use their own files, as long as it is in STL or OBJ format.

The simulation parameters available for the user to configure are:
- `num_iterations`: Set the maximum number of iterations for the iterative algorithm to converge.
- `resolution`: Controls the resolution of the meshing that is done prior to the simulation. The higher the resolution, the finer the meshing. Possibilities: "high", "medium", "low".

Moreover, the hardware and interaction are configured with the usual general parameters - `machine_group`, `run_async`, `n_cores`.
Launching a simulation returns a task object, which can be used to verify the status of the simulation, get the simulation outputs and access post-processing tools. See more in the [Tasks section](inductiva/tasks/README.md).


### Example:

To run this example you need an object of OBJ or STL format. See below how to fetch one that we have available in our Github repository.

Do not forget to insert your API Key (get one by filling this [form](https://docs.google.com/forms/d/e/1FAIpQLSflytIIwzaBE_ZzoRloVm3uTo1OQCH6Cqhw3bhFVnC61s7Wmw/viewform?usp=sf_link)).

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Url to a test an object in Inductiva Github repository
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
    num_iterations=50, resolution="low")

# Download the simulation output to your local machine.
output = task.get_output()
```

This last step is essential to download the simulation files into your local machine 
and apply the post-processing and visualizations to understand the aerodynamics of the object.

## Output and Post-Processing

The simulation output will contain several sub-directories with files that can be post-processed to extract relevant metrics, such as:
- Pressure field: Pressure field of the airflow over the object.
- Streamlines: Streamlines of the airflow over the domain and interacting with the object.
- Flow slice: Slice of the flow that represents the airflow over the domain and interacting with the object.
- Force coefficients: Force coefficients that represent the forces acting on the object. These are the drag and lift coefficients.

### Post-processing Tools

Since at times the output files can be large, the user can choose to get only some of the default outputs that are post-processed
remotely by just doing `output=task.get_output()` (e.g., pressure field, streamlines and flow_slice). Otherwise, the user can choose to get all the outputs with `output=task.get_output(all_files=True)` and post-process as he wishes locally. 

In any case, to obtain these metrics the user can use the following post-processing tools available through the `output` object:
- `get_output_mesh`: Returns two pyvista meshes, one of the airflow over the entire domain and the other of data on the object.
    These meshes contain the results of the simulation, which are breakdown and visualized with the next tools.
- `get_object_pressure_field`: Returns a pyvista mesh with the pressure field of the airflow over the object.
- `get_streamlines`: Input parameters - `max_time`, `n_points`, `initial_step_length`, `source_radius`, `source_center`. Returns a pyvista mesh of the streamlines with pressure and velocity components of airflow over the domain.
- `get_flow_slice`: Input parameters - `plane`, `origin`. Returns a pyvista mesh of the flow slice with pressure and velocity components.
- `get_force_coefficients`: Returns the force coefficients that represent the forces acting on the object. These are the drag and lift coefficients.

#### Visualizations

For the pressure field, streamlines and flow slices there are easy-to-use visualizations, which have some configuration parameters (see the examples below for a general overview of these parameters). To render any of these objects one does `streamlines.render()`. **Remark: ** For those running on Google Colab or a headless server, further extra dependencies are required. Install them with 

```
!apt install libgl1-mesa-glx xvfb
```

Moreover, change the parameter `virtual_display` in the render method to `True`, and pass a path to save the result, like 

```python
streamlines.render(virtual_display=True, save_path="streamlines.png")
```

The next example fetches the default post-processed files and renders the respective visualizations. The last one, fetches the full simulation output, and with the same interface the user can configure the post-processing tools to extract the metrics he desires and render the visualizations.

### Example for default files:

```python
# Get the default files from WindTunnel simulation
output = task.get_output()

# Get a pressure_field mesh
pressure_field = output.get_object_pressure_field()

# Render
pressure_field.render()
```

<img src="/resources/media/openfoam/default_pressure_field.png" width="400" height="300" />

```python
# Get a mesh of the streamlines 
streamlines = output.get_streamlines()

# Render the streamlines
streamlines.render(physical_field="velocity",
                   flow_cmap="viridis",
                   view="isometric",
                   streamline_radius=0.1
                   save_path="default_streamlines.png")
```

<img src="/resources/media/openfoam/default_streamlines.png" width="400" height="300" />


```python
# Get a mesh of the flow slice
flow_slice = output.get_flow_slice(plane="xz")

# Render the pressure field over the slice
flow_slice.render(physical_field="pressure",
                  flow_cmap="viridis",
                  save_path="default_flow_slice.png")
```

<img src="/resources/media/openfoam/default_flow_slice.png" width="400" height="300" />

### Example for general Post-processing

```python
# Get all the WindTunnel simulation files
output = task.get_output(all_files=True)

# Get a pressure field mesh
pressure_field = output.get_object_pressure_field()

# Render
pressure_field.render(save_path="pressure_field.png")
```

<img src="/resources/media/openfoam/pressure_field.png" width="400" height="300" />

```python
# Get the streamlines mesh
streamlines = output.get_streamlines(max_time=200,
                                     n_points=200,
                                     initial_step_length=0.1,
                                     source_radius=0.5,
                                     source_center=[-3, 0, 1])

# Render the streamlines
streamlines.render(physical_field="pressure",
                   flow_cmap="viridis",
                   view="isometric",
                   streamline_radius=0.1,
                   save_path="streamlines.png")
```

<img src="/resources/media/openfoam/streamlines.png" width="400" height="300" />

```python

# Get the flow slice mesh
flow_slice = output.get_flow_slice(plane="xz",
                                   origin=(0,0,0))

# Render the velocity field over the flow slice.
flow_slice.render(physical_field="velocity",
                  flow_cmap="Blues",
                  save_path="flow_slice.png")
```

<img src="/resources/media/openfoam/flow_slice.png" width="400" height="300" />
