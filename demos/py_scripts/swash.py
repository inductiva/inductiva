"""Sample usage of SWASH simulator via the API."""
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
flags.DEFINE_string("sim_config_filename",
                    None,
                    "Name of the input .sws file.",
                    required=True)
flags.DEFINE_string("output_dir", None,
                    "Directory where the outputs will be stored.")

flags.DEFINE_string("working_dir", None,
                    "Directory to which paths are relative to.")


def main(_):
    """Run a SWASH simulation using user-provided input files."""

    inductiva.api_url = FLAGS.api_url
    inductiva.working_dir = FLAGS.working_dir

    swash_sim = inductiva.simulators.SWASH()

    task = swash_sim.run(
        input_dir=FLAGS.sim_dir,
        sim_config_filename=FLAGS.sim_config_filename,
    )

    output = task.get_output(output_dir=FLAGS.output_dir)
    logging.info("Task ID %s", task.id)
    logging.info("Outputs stored in %s", output)


if __name__ == "__main__":
    app.run(main)
