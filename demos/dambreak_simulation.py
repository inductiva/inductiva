"""Sample usage of SPlisHSPlasH simulation via API.
"""
import inductiva

from absl import logging

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.update_config(address="http://192.168.1.50:8000", output_dir="output")

    scenario = inductiva.fluids.DamBreak(fluid=inductiva.fluids.WATER,
                                         fluid_dimensions=[0.05, 0.8, 0.8])
    simulation_output = scenario.simulate()
    simulation_output.render()
