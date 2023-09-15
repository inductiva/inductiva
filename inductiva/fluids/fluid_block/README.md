# Fluid Block

This scenario simulates the motion of a fluid in a cubic tank. 
The fluid starts as a cube and then flows under the effect of gravity.
The simulation is performed using the [Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)
method. Currently, this scenario runs with the following simulators: [SPlisHSPlasH](https://github.com/inductiva/inductiva/blob/main/inductiva/simulators/splishsplash.py) and [DualSPHysics](https://github.com/inductiva/inductiva/blob/main/inductiva/simulators/dualsphysics.py).

This scenario allows users to play with the different properties of fluids, namely, the density and the kinematic viscosity of the fluid.
Further, users can select the size of the initial fluid block, the initial position and the initial velocity. 
After the scenario has been initialized, users are ready to launch the simulation. They can configure the simulation further with the following parameters: the particule radius (float), the simulation time (float), adaptive time step (float), the time step (float) and the output time step (float).

### Example

```python
import inductiva

inductiva.api_key = "ADD_API_KEY_HERE"

# Initialize the scenario
scenario = inductiva.fluids.FluidBlock(
    density=1000,
    kinematic_viscosity=0.001,
    dimensions=[0.4, 0.3, 0.4])

# Default simulator is Dualsphysics
task = scenario.simulate(simulation_time=2)

output = task.get_output()
```

To visualize the results, users will need to install the extra dependencies of the fluid package with `pip install --upgrade "inductiva[fluids_extra]"`. 
Then run the following code:

```python
output.render()
```

<p align="center">
  <img src="https://github.com/inductiva/inductiva/assets/" alt="Centered Image" width="350" height="250">

Remark: For those running on Google Colab or a headless server, further extra dependencies are required. Install them with `!apt install libgl1-mesa-glx xvfb`. Moreover, activate the parameter `virtual_display=True` in the render method, like `output.render(virtual_display=True)`.

Further, ou can alter the `color` of the particles and `fps` of the rendering!
