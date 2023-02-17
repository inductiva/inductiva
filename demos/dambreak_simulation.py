"""Sample usage of SPlisHSPlasH simulation via API."""
import inductiva

from absl import logging
from absl import flags
from absl import app

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")


def main(_):
    inductiva.api_url = FLAGS.api_url

    scenario = inductiva.fluids.DamBreak(fluid=inductiva.fluids.WATER,
                                         fluid_dimensions=[0.05, 0.8, 0.8])
    simulation_output = scenario.simulate()
    simulation_output.render()


if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)
    app.run(main)
