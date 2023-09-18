# Dam Break

This scenario simulates the motion of a fluid in a cubic tank under the effect of gravity.
The simulation is performed using the [Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)
method. Currently, this scenario runs with the following simulators: [SPlisHSPlasH](https://github.com/inductiva/inductiva/blob/main/inductiva/simulators/splishsplash.py) and [DualSPHysics](https://github.com/inductiva/inductiva/blob/main/inductiva/simulators/dualsphysics.py).

In this scenario, users can pick different fluids (like Water, Honey, Olive oil, Beer, Jet fuel, and more, listed [here](https://github.com/inductiva/inductiva/blob/main/inductiva/fluids/fluid_types.py)), adjust the cube's size, and set its initial position.

After the scenario has been initialized, users are ready to launch the simulation. They can configure the simulation further with the following parameters: the resolution - "low", "medium", "high" - and the simulation time.

### Example

Do not forget to insert your API Key (get one by filling this [form](https://docs.google.com/forms/d/e/1FAIpQLSflytIIwzaBE_ZzoRloVm3uTo1OQCH6Cqhw3bhFVnC61s7Wmw/viewform?usp=sf_link)).

```python
import inductiva

inductiva.api_key = "ADD_API_KEY_HERE"

# Initialize the scenario
scenario = inductiva.fluids.DamBreak()

# Default simulator is DualSPHysics
task = scenario.simulate()

output = task.get_output()
```

To visualize the results, users will need to install some extra dependencies of the fluid package with `pip install --upgrade "inductiva[fluids_extra]"`. 
Then run the following code:

```python
output.render()
```

<p align="center">
  <img src="https://github.com/inductiva/inductiva/blob/f52d0a733276996e02fdde942a4974c0a75d5038/resources/media/dam_break.gif" alt="Centered Image" width="400" height="300">

**Remark:** For those running on Google Colab or a headless server, further extra dependencies are required. Install them with `!apt install libgl1-mesa-glx xvfb`. Moreover, change the parameter `virtual_display` in the render method to `True`, like `output.render(virtual_display=True)`.

Further, one can alter the `color` of the particles and `fps` of the rendering!
