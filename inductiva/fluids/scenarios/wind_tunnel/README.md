# Wind Tunnel Scenario

This scenario models the aerodynamics of an object inside a virtual wind tunnel for a given airflow velocity. Air is injected on a side wall of the wind tunnel, the flow changes according to the structure of the object and leaves through an outlet on the other side. The system is modelled with the steady-state equations for incompressible flow and the $k-\epsilon$ turbulence models.

To initialize the scenario, the user can define the following parameters:
- Dimensions of the wind tunnel in meters, e.g. `{"x": [-5, 15], "y": [-5, 5], "z": [0, 8]}`.
- Airflow velocity in m/s, e.g. `[30, 0, 0]`. Notice that the airflow is injected from the negative x-direction.

Now, the user is ready to simulate the steady state. Here, the user chooses the object to be inserted inside the wind tunnel. This object is defined with a geometry file in STL or OBJ format (currently tested). Notice, that appropriate meshing will be done during the simulation.

The user can further define the following simulation parameters:
- simulator: The simulator to be used for the simulation. Currently, only `OpenFOAM` is supported.
- num_iterations: Number of iterations for the steady-state simulation.
- resolution: Resolution of the meshing. The higher the resolution, the finer the meshing.

Moreover, the hardware and interaction are configured with the usual general parameters - `machine_group`, `run_async`, `n_cores`.
Launching a simulation returns a task object, which can be used to verify the status of the simulation, get the simulation outputs and access instantly post-processing tools. See more in [Tasks](inductiva/tasks/README.md).

### Example:

```python
import inductiva

# Initialize the scenario
scenario = inductiva.fluids.scenarios.WindTunnel(
    flow_velocity=[30, 0, 0],
    dimensions={"x": [-5, 15], "y": [-5, 5], "z": [0, 8]})

# Run a simulation
task = scenario.simulate(
    object_path="f1.obj",
    num_iterations=1000, resolution=0.5,
    run_async=True, n_cores=4)

# Get the simulation output on your local machine.
output = task.get_output()
```

## Output and Post-Processing:

The `WindTunnel` scenario allows users to extract metrics about the object under the established flow velocity.
The simulation outputs of the wind tunnel are:
- `pressure_field`: Pressure field of the airflow over the object.
- `streamlines`: Streamlines of the airflow over the domain and interacting with the object.
- `flow_slice`: Slice of the flow that represents the airflow over the domain and interacting with the object.
- `force_coefficients`: Force coefficients that represent the forces acting on the object. These are the drag and lift coefficients.

Since at times the output files can be large, the user can choose to get only some of the default outputs that are post-processed
on run-time by just doing `task.get_output()`. Otherwise, the user can choose to get all the outputs with `task.get_output(full=True)` and post-process as he wishes locally. 

In any case, to obtain the metrics above the user can use the following post-processing tools:
- `get_output_mesh`: Returns the mesh of the airflow over the entire domain and data on the object.
- `get_object_pressure_field`: Returns the pressure field of the airflow over the object.
- `get_streamlines`: Input parameters - `max_time`, `n_points`, `initial_step_length`, `source_radius`, `source_center`. Returns a mesh of the streamlines with pressure and velocity components of airflow over the domain.
- `get_flow_slice`: Input parameters - `plane`, `origin`. Returns a mesh of the flow slice with pressure and velocity components.
- `get_force_coefficients`: Returns the force coefficients that represent the forces acting on the object. These are the drag and lift coefficients.

For the pressure field, streamlines and flow slices there are easy-to-use visualizations. See the example below.

### Example:

```python

# Get the full output from WindTunnel simulation
output = task.get_output(full=True)

pressure_field = output.get_object_pressure_field()
pressure_field.render_frame()
```

```python
# Get the streamlines and visualize them
streamlines = output.get_streamlines(max_time=200,
                                     n_points=200,
                                     initial_step_length=0.1,
                                     source_radius=0.5,
                                     source_center=[-3, 0, 1])
streamlines.render_frame()
```

```python

flow_slice = output.get_flow_slice(plane="xz", origin=[0, 0, 0])
flow_slice.render_frame()
```


