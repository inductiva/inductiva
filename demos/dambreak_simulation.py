"""Sample usage of SPlisHSPlasH simulation via API.
"""
import inductiva

from absl import logging
from absl import flags

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")

if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.update_config(FLAGS.api_url)

    scenario = inductiva.fluids.DamBreak(fluid=inductiva.fluids.WATER,
                                         fluid_dimensions=[0.05, 0.8, 0.8])
    simulation_output = scenario.simulate()
    simulation_output.render()
