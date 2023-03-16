"""Sample usage of dam break simulation via API."""
import time

from absl import logging
from absl import flags
from absl import app

import inductiva
from inductiva.fluids.scenarios import DualSPHysicsParameters
import inductiva_utils

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")
flags.DEFINE_list("fluid_dimensions", [0.2, 1, 1],
                  "Dimensions of the fluid column.")
flags.DEFINE_list("fluid_position", [0.0, 0.0, 0.0],
                  "Position of the fluid column in the tank.")
flags.DEFINE_enum("resolution", "medium", ["high", "medium", "low"],
                  "Sets the fluid resolution to simulate.")
flags.DEFINE_enum("engine", "DualSPHysics", ["DualSPHysics", "SPlisHSPlasH"],
                  "Sets the fluid resolution to simulate.")
flags.DEFINE_string("output_dir", None,
                    "Destination directory for output files.")
flags.DEFINE_float("simulation_time", 1, "Simulation time in seconds.")
flags.DEFINE_string("device", "cpu",
                    "Device in which device the simulation will run.")


def main(_):
    """Run a dam break simulation via the API."""
    inductiva.api_url = FLAGS.api_url

    time_start = time.perf_counter()

    scenario = inductiva.fluids.DamBreak(
        fluid=inductiva.fluids.WATER,
        fluid_dimensions=inductiva_utils.flags.cast_list_to_float(
            FLAGS.fluid_dimensions),
        fluid_position=inductiva_utils.flags.cast_list_to_float(
            FLAGS.fluid_position))

    _ = scenario.simulate(output_dir=FLAGS.output_dir,
                          resolution=FLAGS.resolution,
                          engine=FLAGS.engine,
                          simulation_time=FLAGS.simulation_time,
                          device=FLAGS.device)

    # Note: video rendering only works with SPlisHSPlasH for now
    # simulation_output.render()

    logging.info("Local time: %s", time.perf_counter() - time_start)


if __name__ == "__main__":
    app.run(main)
