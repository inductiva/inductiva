"""Demo simulation of a fluid tank via API."""
import time

from absl import logging
from absl import flags
from absl import app

import inductiva
from inductiva.fluids.scenarios.fluid_tank import FluidTank
from inductiva.fluids.shapes import Cylinder

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")
flags.DEFINE_float("tank_radius", 0.5, "Tank radius, in meters.")
flags.DEFINE_float("tank_height", 1, "Tank height, in meters.")
flags.DEFINE_float("fluid_level", 0.5,
                   "Initial fluid level inside the tank, in meters.")

flags.DEFINE_string("output_dir", "output",
                    "Destination directory for output files.")
flags.DEFINE_enum("device", "cpu", ["cpu", "gpu"],
                  "Device in which device the simulation will run.")


def main(_):
    """Run a fluid tank simulation via API."""
    inductiva.api_url = FLAGS.api_url

    time_start = time.perf_counter()

    scenario = FluidTank(
        fluid=inductiva.fluids.WATER,
        shape=Cylinder(radius=FLAGS.tank_radius, height=FLAGS.tank_height),
        fluid_level=FLAGS.fluid_level,
    )

    scenario.simulate(simulator=inductiva.simulators.SPlisHSPlasH(),
                      device=FLAGS.device)

    logging.info("Local time: %s", time.perf_counter() - time_start)


if __name__ == "__main__":
    app.run(main)
