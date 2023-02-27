"""Sample usage of SPlisHSPlasH simulation via the API."""
from absl import logging
from absl import flags
from absl import app

import inductiva

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://localhost:8000",
                    "Base URL of the Inductiva API.")
flags.DEFINE_string("sim_dir", "testcases/a11stwav",
                    "Directory with the simulation inputs.")
flags.DEFINE_string("input", "a11stw01", "Name of the input .sws file.")
flags.DEFINE_string("output_dir", "testing",
                    "Directory where the outputs will be stored.")
flags.DEFINE_integer("n_cores", 4, "Number of cores to use.")


def main(_):
    inductiva.api_url = FLAGS.api_url

    inductiva.swash.run_simulation(
        FLAGS.sim_dir,
        n_cores=FLAGS.n_cores,
        input=FLAGS.input,
        output_dir=FLAGS.output_dir,
    )


if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)
    app.run(main)
