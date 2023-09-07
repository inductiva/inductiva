### Fluid tank

This scenario simulates the motion of a fluid in a cubic or cylindrical tank.
Fluid is injected into the tank via an inlet located at the top of the tank and
flows out of the tank via an outlet located at the bottom of the tank. The
motion of the fluid is controlled by gravity. The simulation is performed using
the [Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)
method.

#### Example

Initialize the scenario:

```python
from inductiva import fluids
scenario = fluids.scenarios.FluidTank(shape=fluids.shapes.Cylinder(),
                                      fluid=fluids.WATER,
                                      fluid_level=0.5)
```

The user can specify the fluid (e.g. water, honey or oil), the fluid level (in
meters), as well as the tank shape and its inlet and outlet properties.

Run the simulation:

```python
task = scenario.simulate(simulation_time=5,
                           output_time_step=0.1,
                           resolution="medium")

output = task.get_output()
```

The user can specify the total simulation time and the time step between outputs
(all in seconds). The user can also specify the resolution of the simulation
(low, medium or high).

Visualize the results:

```python
output.render()
```

![Fluid tank simulation.](resources/media/fluid_tank.gif)