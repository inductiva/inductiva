import os

from inductiva.types import DirPath
from IPython.display import HTML
from inductiva_data.data import ParticleDataReader
from inductiva_data import visualizers


class SimulationOutput:

    def __init__(self, sim_output_path: DirPath) -> None:
        self.sim_output_dir=sim_output_path.path

    def render(self):
        # Read simulation particle data
        reader = ParticleDataReader()

        particle_data = reader.read_dir(os.path.join(self.sim_output_dir, "hdf5"))

        visaulizer = visualizers.TimeVaryingParticleData3DScatterVisualizer(
            data=particle_data,
            x_quantity="x",
            y_quantity="y",
            z_quantity="z",
        )
        movie_dir = os.path.join(self.sim_output_dir, "movie.mp4")
        visaulizer.create_time_movie(movie_dir)

        return HTML("""
            <video alt="test" controls>
                <source src="movie_dir" type="video/mp4">
            </video>
        """)


