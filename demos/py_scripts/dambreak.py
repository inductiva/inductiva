"""Sample usage of dam break simulation via API."""
import time

from absl import logging
from absl import flags
from absl import app

import inductiva
from inductiva.simulators import SplishSplash, DualSPHysics
from inductiva.fluids.fluid_types import WATER
# from inductiva.fluids.fluid_types import get_fluid_color
from inductiva.utils.flags import cast_list_to_float

SIMULATORS = {
    "splishsplash": SplishSplash,
    "dualsphysics": DualSPHysics,
}

FLAGS = flags.FLAGS

flags.DEFINE_string("api_url", "http://api.inductiva.ai",
                    "Base URL of the Inductiva API.")
flags.DEFINE_list("fluid_dimensions", [0.2, 1, 1],
                  "Dimensions of the fluid column.")
flags.DEFINE_list("fluid_position", [0.0, 0.0, 0.0],
                  "Position of the fluid column in the tank.")
flags.DEFINE_enum("resolution", "low", ["high", "medium", "low"],
                  "Sets the fluid resolution to simulate.")
flags.DEFINE_enum("simulator", "splishsplash", SIMULATORS.keys(),
                  "Fluid simulator to use.")
flags.DEFINE_string("output_dir", None,
                    "Destination directory for output files.")
flags.DEFINE_float("simulation_time", 1, "Simulation time in seconds.")
flags.DEFINE_string("device", "gpu",
                    "Device in which device the simulation will run.")


def main(_):
    """Run a dam break simulation via the API."""
    inductiva.api_url = FLAGS.api_url

    time_start = time.perf_counter()

    scenario = inductiva.fluids.DamBreak(
        fluid=WATER,
        dimensions=cast_list_to_float(FLAGS.fluid_dimensions),
        position=cast_list_to_float(FLAGS.fluid_position))

    simulator_cls = SIMULATORS[FLAGS.simulator]

    task = scenario.simulate(
        simulator=simulator_cls(),
        resolution=FLAGS.resolution,
        simulation_time=FLAGS.simulation_time,
        device=FLAGS.device,
    )

    # Note: video rendering only works with SPlisHSPlasH for now
    # simulation_output.render(color=get_fluid_color(WATER), alpha=0.8)
    _ = task.get_output()

    logging.info("Local time: %s", time.perf_counter() - time_start)


if __name__ == "__main__":
    app.run(main)
