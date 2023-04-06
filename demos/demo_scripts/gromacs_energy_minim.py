"""Sample usage of Gromacs energy minimization via the API."""

from absl import logging
from absl import flags
from absl import app

import inductiva

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://localhost:6969",
                    "Base URL of the Inductiva API.")

flags.DEFINE_string("sim_dir",
                    None,
                    "Directory with the simulation inputs.",
                    required=True)

flags.DEFINE_string("sim_config_filename",
                    None,
                    "Name of the input file.",
                    required=True)

flags.DEFINE_string("protein_filename",
                    None,
                    "Name of the protein file.",
                    required=True)

flags.DEFINE_string("topology_filename",
                    None,
                    "Name of the topology file.",
                    required=True)

flags.DEFINE_string("output_dir",
                    None,
                    "Directory where the outputs will be stored.",
                    required=True)

def main(_):
    """Run a Gromacs energy minimization using user-provided input files."""

    inductiva.api_key = "user-key"
    inductiva.api_url = FLAGS.api_url

    gromacs_sim = inductiva.md.simulators.GROMACS()

    output_path = gromacs_sim.run(input_dir=FLAGS.sim_dir,
                                  sim_config_filename=FLAGS.sim_config_filename,
                                  protein_filename=FLAGS.protein_filename,
                                  topology_filename=FLAGS.topology_filename,
                                  output_dir=FLAGS.output_dir)

    logging.info("Outputs stored in %s", output_path)

if __name__ == "__main__":
    app.run(main)

