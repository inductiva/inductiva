"""Sample usage of CoastalArea scenario via the API."""

from absl import logging
from absl import flags
from absl import app

import inductiva

FLAGS = flags.FLAGS

import inductiva

def main(_):
    bathymetry = inductiva.coastal.Bathymetry.from_random_depths(
        x_range=(0, 100),
        y_range=(0, 100),
        x_num=100,
        y_num=100,
        max_depth=1,
        random_seed=12)

    scenario = inductiva.coastal.CoastalArea(bathymetry=bathymetry,
                                         wave_source_location="W",
                                         wave_amplitude=0.1,
                                         wave_period=5.5)

    task = scenario.simulate(simulation_time=360, output_time_step=1)

    output = task.get_output()

    output.render(movie_path="movie_path.mp4")

if __name__ == "__main__":
    app.run(main)
