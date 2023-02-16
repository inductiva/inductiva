"""Post process SPlisHSPlasH simulation outputs."""
import os

from IPython.display import HTML
from base64 import b64encode

from inductiva.types import DirPath
from inductiva_data.data import ParticleDataReader
from inductiva_data import visualizers


class SimulationOutput:
    """Post process SPlisHSPlasH simulation outputs."""

    def __init__(self, sim_output_path: DirPath) -> None:
        """Initializes a `SimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files."""
        self.sim_output_dir = sim_output_path.path

    def render(self):
        """Generate a simulation movie."""

        # Read simulation particle data
        reader = ParticleDataReader()
        particle_data = reader.read_dir(
            os.path.join(self.sim_output_dir, "hdf5"))

        visualizer = visualizers.TimeVaryingParticleData3DScatterVisualizer(
            data=particle_data,
            x_quantity="x",
            y_quantity="y",
            z_quantity="z",
        )
        movie_path = os.path.join(self.sim_output_dir, "movie.mp4")
        visualizer.create_time_movie(movie_path)

        with open(movie_path, "r") as fp:
            mp4 = fp.read()
        movie_url = "data:video/mp4;base64," + b64encode(mp4).decode()

        return HTML(f"""
            <video alt="test" controls>
                <source src="{movie_url}" type="video/mp4">
            </video>
        """)
