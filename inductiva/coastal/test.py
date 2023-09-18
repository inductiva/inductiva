"""End to end CoastalArea scenario simulation"""

import inductiva
import matplotlib.pyplot as plt

bathymetry = inductiva.coastal.Bathymetry.from_random_depths(x_range=(0, 100),
                                                             y_range=(0, 100),
                                                             x_num=100,
                                                             y_num=100,
                                                             max_depth=50,
                                                             random_seed=12)

bathymetry.plot()
plt.show()

scenario = inductiva.coastal.CoastalArea(bathymetry=bathymetry,
                                         wave_source_location="W",
                                         wave_amplitude=5,
                                         wave_period=5.5)

task = scenario.simulate(simulation_time=30, output_time_step=1)

output = task.get_output()

output.render(movie_path="movie_path.mp4")
