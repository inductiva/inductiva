"""Sample usage of SPlisHSPlasH simulation via the API."""
from absl import logging
from absl import flags
from absl import app

import inductiva

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")
flags.DEFINE_string("sim_dir",
                    None,
                    "Directory with the simulation inputs.",
                    required=True)
flags.DEFINE_string("output_dir",
                    None,
                    "Directory where the outputs will be stored.",
                    required=True)


def main(_):
    inductiva.api_url = FLAGS.api_url

    inductiva.sph.run_simulation(
        FLAGS.sim_dir,
        output_dir=FLAGS.output_dir,
    )


if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)
    app.run(main)
