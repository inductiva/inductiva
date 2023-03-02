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
flags.DEFINE_string("input_filename",
                    None,
                    "Name of the input .sws file.",
                    required=True)
flags.DEFINE_string("output_dir", None,
                    "Directory where the outputs will be stored.")
flags.DEFINE_integer("n_cores", 1, "Number of cores to use.")


def main(_):
    """Run a SPlisHSPlasH simulation using user-provided input files."""

    inductiva.api_url = FLAGS.api_url

    swash_sim = inductiva.fluids.Swash(
        sim_dir=FLAGS.sim_dir,
        input_filename=FLAGS.input_filename,
    )

    output_path = swash_sim.simulate(
        output_dir=FLAGS.output_dir,
        n_cores=FLAGS.n_cores,
    )

    logging.info("Outputs stored in %s", output_path)


if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)
    app.run(main)
