# Wind Tunnel Scenario

This scenario models the aerodynamics of an object inside a virtual wind tunnel for a given airflow velocity. To model the air as an incompressible fluid and simplify the equations of motion under study, we restrict ourselves here to maximum speeds available is 100 m/s. Above this, the model is inaccurate since it does not account for the compressibility of the air at those speeds.
Air is injected through one side of the wind tunnel and exits through an outlet on the opposite side. The airflow within the tunnel is modified according to the structure of the object. The system is modelled with the steady-state equations for incompressible flow and the $k-\epsilon$ turbulence models. The simulation is currently performed with the OpenFOAM simulator.

### Example:

**Note:** To access the visualization methods the user needs to install the extra dependencies for the fluids package with `pip install --upgrade "inductiva[fluids_extra]"`.

This example requires an object of OBJ or STL format, we make one available through the code.

Do not forget to insert your API Key (check the [main page](https://github.com/inductiva/inductiva/tree/main#api-access-tokens) to see how to get one).

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Url to a test an object in Inductiva Github repository
vehicle_url = "https://raw.githubusercontent.com/inductiva/inductiva/main" \
              "/assets/vehicle.obj"
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

# Get a pressure_field mesh
pressure_field = output.get_object_pressure_field()

# Render
pressure_field.render()
```
<img src="/resources/media/openfoam/default_pressure_field.png" width="400" height="300" />

## Scenario configuration

To initialize the scenario, the user can define the following parameters:
- `domain`: Dimensions of the wind tunnel in meters, e.g. `{"x": [-5, 15], "y": [-5, 5], "z": [0, 8]}`.
- `flow_velocity`: Airflow velocity in m/s, e.g. `[30, 0, 0]`. Notice that the airflow is injected from the negative x-direction.

Now, the user is ready to simulate the steady state. Here, the user chooses the object to be inserted inside the wind tunnel. This object is defined with a geometry file in [STL](https://en.wikipedia.org/wiki/STL_(file_format)) or [OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) format. Notice that:
- the object is not re-scaled automatically, so the user must either change the dimensions of the domain or scale the object himself.
- the required meshing step is automatically done on the backend before the simulation starts with [snappyHexMesh tool](https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.4-mesh-generation-with-the-snappyhexmesh-utility) of OpenFOAM.

The simulation parameters available for the user to configure are:
- `num_iterations`: Set the maximum number of iterations for the iterative algorithm to converge.
- `resolution`: Controls the resolution of the meshing that is done prior to the simulation. The higher the resolution, the finer the meshing. Available options: "high", "medium", "low".

Moreover, the hardware is configured through the creation of a `machine_group` to run your simulations - see more in the [Machine Groups section](inductiva/resources/README.md).
Launching a simulation returns a task object, which can be used to verify the status of the simulation, get the simulation outputs and access post-processing tools. See more in the [Tasks section](inductiva/tasks/README.md).

## Output and Post-Processing

The simulation output will contain several sub-directories with files that can be post-processed to extract relevant metrics, such as:
- Pressure field: Pressure field of the airflow over the object.
- Streamlines: Streamlines of the airflow over the domain and interacting with the object.
- Flow slice: Slice of the flow that represents the airflow over the domain and interacting with the object.
- Force coefficients: Force coefficients that represent the forces acting on the object. These are the drag and lift coefficients.

### Post-processing methods

Since at times the simulation output files can be large, e.g. 1 Gb, the user can choose to download only some of the default outputs that are post-processed
remotely by just doing `output=task.get_output()` (e.g., pressure field, streamlines and flow_slice). Alternatively, the user can choose to download all simulation outputs with `output=task.get_output(all_files=True)` and post-process as wishes locally. 

In any case, to obtain these metrics the user can use the following post-processing methods available through the `output` object:
- `get_output_mesh`: Returns two pyvista meshes, one of the airflow over the entire domain and the other of data on the object.
    These meshes contain the results of the simulation, which are breakdown and visualized with the next methods.
- `get_object_pressure_field`: Returns a MeshData object - object with scalar measurements at each point and face of the mesh, as well as rendering capabilities -  with the pressure field over the object mesh. This object can be rendered with `render()`.
- `get_streamlines`: Input parameters - `max_time`, `n_points`, `initial_step_length`, `source_radius`, `source_center`. Returns a `Streamlines` - an object with the mesh of the streamlines and a `render()` method with the vehicle inside the WindTunnel. The streamlines have pressure and velocity components of airflow over the domain.
- `get_flow_slice`: Input parameters - `plane`, `origin`. Returns a `FlowSlice` - an object with the mesh of a slice of the domain with flow properties over the mesh, and a `render()` method with the vehicle inside the WindTunnel. The flow slice has pressure and velocity components of the airflow over the domain.
- `get_force_coefficients`: Returns a list with the force coefficients of the steady state - Drag, Lift, Front Lift and Rear Lift - that represent the forces acting on the object. 

Each visualization method has extra configuration parameters that can be passed to the `render()` method, see in the examples below for an overview for each of them.

**Remark: **
- For those running on Google Colab or a headless server, further extra dependencies are required. Install them with `!apt install libgl1-mesa-glx xvfb`

Moreover, change the parameter `virtual_display` in the render method to `True`, and pass a path to save the result, like 

```python
streamlines.render(virtual_display=True, save_path="streamlines.png")
```

- To save your visualizations, please use the 'q' key to close the plotter. Some operating systems, particularly Windows, may encounter issues when attempting to save a screenshot if you use the exit button in the GUI.

The next example fetches the default post-processed files and renders the respective visualizations. The last one, fetches the full simulation outputs, and with the same interface the user can configure the post-processing methods to extract the metrics he desires and render the visualizations.

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
