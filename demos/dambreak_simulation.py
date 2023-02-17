"""Sample usage of SPlisHSPlasH simulation via API.
"""
from absl import logging
from absl import flags

import inductiva
import inductiva_utils

FLAGS = flags.FLAGS

flags.DEFINE_list("fluid_dimensions", [0.2, 0.8, 0.8], 
                  "Dimensions of the fluid column.")
flags.DEFINE_list("fluid_position", [0.0, 0.0, 0.0], 
                  "Position of the fluid column in the tank.")
flags.DEFINE_float("particle_radius", 0.02, 
                  "Radius of the discretization particles, in meters.")


if __name__ == "__main__":
    logging.set_verbosity(logging.DEBUG)

    inductiva.init(address="http://192.168.1.50:8000", output_dir="output")

    scenario = inductiva.fluids.DamBreak(
        fluid=inductiva.fluids.WATER,
        fluid_dimensions=inductiva_utils.flags.cast_list_to_float(FLAGS.fluid_dimensions),
        fluid_position=inductiva_utils.flags.cast_list_to_float(FLAGS.fluid_position), 
        particle_radius=FLAGS.particle_radius)

    simulation_output = scenario.simulate()
    simulation_output.render(color_quantity="y")
