"""Sample usage of Gromacs energy minimization via the API."""

from absl import logging
from absl import flags
from absl import app

import inductiva

flags.DEFINE_string("input_dir",
                    "demos/config_files/gromacs",
                    "Directory containing the input files.",)
                    # required=True)

flags.DEFINE_string("protein_pdb",
                    "alpline.pdb",
                    "Path to the protein PDB file.")

flags.DEFINE_string("output_dir",
                    None,
                    "Directory where the outputs will be stored.",)
                    # required=True)

FLAGS = flags.FLAGS

def main(_):
    """Run a Gromacs command."""
    inductiva.api_key = "735d2e8d919a5358ad2cd478d5abd9805482cc4b087dac4a7f20a34993f6e411"

    gromacs = inductiva.molecules.simulators.GROMACS()

    example_command = {
        "cmd":
            f"gmx pdb2gmx -f {FLAGS.protein_pdb} -o protein.gro -water tip3p -ff amber99sb-ildn",
        "promps": []
    }

    gromacs.run(input_dir=FLAGS.input_dir,
                commands=example_command,
                output_dir=FLAGS.output_dir)


if __name__ == "__main__":
    app.run(main)
