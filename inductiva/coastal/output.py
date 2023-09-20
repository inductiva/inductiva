"""Post-processing and visualization utilities of the coastal area scenario."""

import os
import tempfile
from typing import Dict, Literal, Optional, Sequence, Tuple
from absl import logging

import numpy as np
try:
    import matplotlib.pyplot as plt
    import scipy
except ImportError:
    plt = None
    scipy = None

from inductiva.types import Path
from inductiva.utils import visualization, optional_deps

# Labels used in SWASH simulation outputs.
QUANTITY_SWASH_LABELS = {
    "water_level": "WATLEV",
    "velocity_x": "VKSI",
    "velocity_y": "VETA",
    "velocity_magnitude": "VMAG",
}

# Labels used in plots.
QUANTITY_PLOT_LABELS = {
    "water_level": "Water level",
    "velocity_x": "Velocity x-component",
    "velocity_y": "Velocity y-component",
    "velocity_magnitude": "Velocity magnitude",
}

# Units used in plots.
QUANTITY_UNITS = {
    "water_level": "m",
    "velocity_x": "m/s",
    "velocity_y": "m/s",
    "velocity_magnitude": "m/s",
}

GRID_POSITIONS_FILE_NAME = "grid_positions.mat"


@optional_deps.needs_coastal_extra_deps
def read_swash_output_file(file_path: str) -> Dict[str, np.ndarray]:
    """Reads a SWASH output file.

    The output is read as a dictionary, where the keys are the names of the
    variables and the values are the corresponding arrays.
    """

    try:
        # Read file contents
        data_dict = scipy.io.loadmat(file_path)

    except:
        raise ValueError(f"Error reading file '{file_path}'.") from None

    return data_dict


@optional_deps.needs_coastal_extra_deps
def read_swash_grid_positions_file(file_path: str) -> Tuple[np.ndarray]:
    """Reads a SWASH output file containing grid positions."""

    # Read grid positions
    grid_positions_dict = read_swash_output_file(file_path)

    # Get grids with x and y coordinates
    x_grid = grid_positions_dict["Xp"].transpose()
    y_grid = grid_positions_dict["Yp"].transpose()

    # Get arrays with x and y coordinates
    x_array = x_grid[:, 0]
    y_array = y_grid[0, :]

    return x_array, y_array


@optional_deps.needs_coastal_extra_deps
def read_swash_grid_quantity_file(file_path: str, quantity: str):
    """Reads a SWASH output file containing a quantity defined over a grid.

    SWASH outputs are stored in .mat files. Here, we first read the contents
    of one of these files as a dictionary, and then convert the dictionary
    to a list of times and a list of data arrays.

    # The list of times is built from the keys of the dictionary:
    # the keys end with a number, in format 123456_789, identifying the time
    # instant 12 hours, 34 minutes, 56 seconds and 789 milliseconds; we
    # strip the times from the keys and convert them to seconds.

    # The list of data arrays is built from the values of the dictionary:
    # the values correspond to the data values at each of the time instants
    # represented by the corresponding keys; we read each value as a 2D
    # tensor with indices (x, y).
    """

    # Read data file
    data_dict = read_swash_output_file(file_path)

    # Convert data dictionary to lists of times and data arrays
    time_list = []
    data_list = []

    for key, data_values in data_dict.items():
        # Check if key starts with the quantity label attributed by SWASH
        if key.startswith(QUANTITY_SWASH_LABELS[quantity].capitalize()):
            # Parse key string to get time in milliseconds
            time_ms = key[len(QUANTITY_SWASH_LABELS[quantity]):].replace(
                "_", "")

            # Convert time to seconds
            time_s = float(time_ms[:2]) * 3600. + \
                     float(time_ms[2:4]) * 60. + \
                     float(time_ms[4:6]) + \
                     float(time_ms[6:]) / 1000.

            # The values are transposed in the file, so we transpose them
            # again for consistency
            transposed_data_values = data_values.transpose()

            time_list.append(time_s)
            data_list.append(transposed_data_values)

    return time_list, data_list


@optional_deps.needs_coastal_extra_deps
def _render_quantity_grid_data(x_array: np.ndarray,
                               y_array: np.ndarray,
                               time_list: Sequence[float],
                               data_list: Sequence[np.ndarray],
                               quantity: str,
                               movie_path: str,
                               fps: int = 10,
                               cmap: Optional[str] = None,
                               clim: Optional[Sequence[float]] = None):
    """Renders temporal evolution of a quantity defined over a grid."""

    with tempfile.TemporaryDirectory() as tmp_dir:

        axes_limits = [
            x_array.min(),
            x_array.max(),
            y_array.min(),
            y_array.max(),
        ]

        if clim is None:
            clim = [
                np.min([np.min(data) for data in data_list]),
                np.max([np.max(data) for data in data_list]),
            ]

        for idx, (time, data) in enumerate(zip(time_list, data_list)):
            fig = plt.figure()
            ax = fig.add_subplot()

            im = ax.imshow(
                np.transpose(data),
                origin="lower",
                cmap=cmap,
                clim=clim,
                extent=axes_limits,
                aspect="equal",
            )

            plt.colorbar(
                im,
                ax=ax,
                label=(f"{QUANTITY_PLOT_LABELS[quantity]} "
                       f"[{QUANTITY_UNITS[quantity]}]"),
            )

            ax.set(
                xlabel="x [m]",
                ylabel="y [m]",
                title=f"Time: {time:.2f} [s]",
            )

            fig.savefig(os.path.join(tmp_dir, f"{idx:06d}.png"))
            plt.close(fig)

        visualization.create_movie_from_frames(frames_dir=tmp_dir,
                                               movie_path=movie_path,
                                               fps=fps)


class CoastalAreaOutput:
    """Heat sink simulation output."""

    def __init__(self, sim_output_path: Path):
        """Initializes a `HeatSinkOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
        """

        self.sim_output_path = sim_output_path
        self.check_stability()

    @optional_deps.needs_coastal_extra_deps
    def check_stability(self):

        _, data_list = read_swash_grid_quantity_file(
            file_path=os.path.join(self.sim_output_path, "water_level.mat"),
            quantity="water_level",
        )

        # Check if water level is stable
        inital_max_water_level = np.max(data_list[0])
        final_max_water_level = np.max(data_list[-1])

        if final_max_water_level > 2 * inital_max_water_level:
            logging.info("Simulation can show unstable results.\n"
                         "Check the simulation visualization to corroborate.\n"
                         "If the results are unstable, try increasing "
                         "the simulation resolution.")

    @optional_deps.needs_coastal_extra_deps
    def render(
        self,
        quantity: Literal[
            "water_level",
            "velocity_x",
            "velocity_y",
            "velocity_magnitude",
        ] = "water_level",
        movie_path: Path = "movie.mp4",
        fps: int = 5,
        cmap: str = "viridis",
        clim: Optional[Sequence[float]] = None,
    ):
        """Renders temporal evolution of a physical quantity.

        A movie is produced representing the temporal evolution of a physical
        quantity exported in simulations of the coastal area scenario. Available
        quantities are: water_level, velocity_x, velocity_y, velocity_magnitude.

        The movie is saved in the `movie_path` location.

        Args:
            quantity: Quantity to render. Available quantities are: water_level,
              velocity_x, velocity_y, velocity_magnitude.
            movie_path: Path to the movie file.
            fps: Number of frames per second to use in the movie.
            cmap: Colormap used to represent the quantity.
            clim: Colorbar limits. When not specified, the limits are
              automatically determined from the data.
        """

        x_array, y_array = read_swash_grid_positions_file(
            os.path.join(self.sim_output_path, GRID_POSITIONS_FILE_NAME))

        time_list, data_list = read_swash_grid_quantity_file(
            file_path=os.path.join(self.sim_output_path, quantity + ".mat"),
            quantity=quantity,
        )

        logging.info("Starting to render video from simulation data...")

        _render_quantity_grid_data(
            x_array=x_array,
            y_array=y_array,
            time_list=time_list,
            data_list=data_list,
            quantity=quantity,
            movie_path=movie_path,
            fps=fps,
            cmap=cmap,
            clim=clim,
        )

        logging.info("Writing mp4 file to %s.", movie_path)
