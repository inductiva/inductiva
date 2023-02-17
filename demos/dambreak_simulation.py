"""Sample usage of SPlisHSPlasH simulation via API.
"""
import inductiva

from absl import logging

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.init(address="http://192.168.1.50:8000", output_dir="output")

    scenario = inductiva.fluids.DamBreak(fluid=inductiva.fluids.WATER,
                                         fluid_dimensions=[0.2, 0.8, 0.8],
                                         fluid_position=[0, 0.0, 0.4], particle_radius=0.01)
    simulation_output = scenario.simulate()
    simulation_output.render(color_quantity="y")
