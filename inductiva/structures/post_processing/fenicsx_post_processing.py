"""Post-process DeformablePlate FEniCSx simulation outputs."""

import os

from typing import List, Optional

try:
    import pyvista as pv
except ImportError:
    pv = None

from inductiva import structures, types, utils

SOLVER_OUTPUT_FILENAME = "results.xdmf"
VALID_FIELD_NAMES = [
    "displacement_x", "displacement_y", "displacement_z", "stress_xx",
    "stress_xy", "stress_xz", "stress_yx", "stress_yy", "stress_yz",
    "stress_zx", "stress_zy", "stress_zz", "strain_xx", "strain_xy",
    "strain_xz", "strain_yx", "strain_yy", "strain_yz", "strain_zx",
    "strain_zy", "strain_zz"
]


class DeformablePlateOutput:
    """Post-process DeformablePlate FEniCSx simulation outputs."""

    def __init__(self, sim_output_path: types.Path):
        """Initializes a 'DeformablePlateOutput' object.

        Args:
            sim_output_path: Path to simulation output files.
        """
        self.sim_output_path = sim_output_path

        # Define the XDMF file path for the solver output file
        file_path = os.path.join(self.sim_output_path, SOLVER_OUTPUT_FILENAME)

        # Remove time values from the XDMF file
        structures.utils.remove_time_values_from_xdmf(file_path, file_path)

    @utils.optional_deps.needs_structures_extra_deps
    def render(self,
               field_names: Optional[List[str]] = None,
               show_edges: bool = True,
               colormap: str = "jet",
               off_screen: bool = True,
               scalar_bar: bool = True,
               background_color: str = "white",
               transparent_background: bool = True,
               virtual_display: bool = True) -> None:
        """Render the simulation fields as images.

        Args:
            field_names (list, optional): List of field names to render. 
              Default is None. The available options for elements in the list
              include "displacement_x," "displacement_y," "displacement_z," 
              "stress_xx", "stress_xy", "stress_xz", "stress_yx", "stress_yy",
              "stress_yz", "stress_zx", "stress_zy", "stress_zz", "strain_xx",
              "strain_xy", "strain_xz", "strain_yx", "strain_yy", "strain_yz",
              "strain_zx", "strain_zy", "strain_zz" and "von_mises".
            show_edges (bool, optional): Whether to display mesh edges. Default
              is True.
            colormap (str, optional): The colormap to apply to the field data.
              Default is "jet".
            off_screen (bool, optional): Whether to render the visualization 
              off-screen. Set to True for non-interactive rendering. Default is 
              True.
            scalar_bar (bool, optional): Whether to include a scalar bar legend
              in the visualization. Default is True.
            background_color (str, optional): The background color of the
              visualization. Default is "white".
            transparent_background (bool, optional): Whether the background of
              the saved images should be transparent. Default is True.
            virtual_display: Whether to use a virtual display to render the
              field.
        """
        # Define the XDMF file path for the solver output file
        file_path = os.path.join(self.sim_output_path, SOLVER_OUTPUT_FILENAME)

        # Read simulation output data
        reader = pv.get_reader(file_path)
        pv_multiblock = reader.read()

        # Get the list of available keys (field names) in the PyVista dataset
        keys = list(pv_multiblock.keys())

        # Select the PyVista dataset representing the mesh
        pv_dataset = pv_multiblock["mesh"]

        if field_names is None:
            field_names = ["von_mises"]

        # Iterate through the specified field names to render
        for field_name in field_names:

            if field_name != "von_mises" and field_name in VALID_FIELD_NAMES:
                # Split the field_name into the main field and subfield
                # (e.g., stress_xx -> field='stress', subfield='xx')
                field, subfield = field_name.split("_")

                # Find the ID of the main field in the dataset keys
                field_id = keys.index(field)

                # Map the subfield string to an ID
                # (e.g., 'xx' -> 0, 'xy' -> 1, 'xz' -> 2)
                subfield_id = {
                    "x": 0,
                    "y": 1,
                    "z": 2,
                    "xx": 0,
                    "xy": 1,
                    "xz": 2,
                    "yx": 3,
                    "yy": 4,
                    "yz": 5,
                    "zx": 6,
                    "zy": 7,
                    "zz": 8
                }[subfield]

                # Extract the specific subfield from the dataset
                # Get the field key from the list of keys
                field_key = keys[field_id]
                # Access the multiblock data using the field key
                multiblock_data = pv_multiblock[field_key]
                # Get the subfield array
                field_array = multiblock_data.get_array(0)[:, subfield_id]

            elif field_name == "von_mises":
                # For the special case of "von_mises" field, find its index in
                # the dataset keys
                field_id = keys.index(field_name)
                # Get the field key from the list of keys
                field_key = keys[field_id]
                # Access the multiblock data using the field key
                multiblock_data = pv_multiblock[field_key]
                # Extract the "von_mises" field array from the dataset
                field_array = multiblock_data.get_array(0)

            else:
                valid_field_names_str = ", ".join(VALID_FIELD_NAMES)
                raise ValueError(
                    f"Invalid field_name {field_name}. "
                    f"Valid field names are: {valid_field_names_str}.")

            # Set the field_array in the point_data of the dataset
            pv_dataset.point_data[field_name] = field_array

        # Create the directory to save the field visualization images
        os.makedirs(os.path.join(self.sim_output_path, "fields"))

        # Define the path to save the field visualization as an image
        field_dir = os.path.join(self.sim_output_path, "fields")

        # Create a 2D fields visualizations from the PyVista dataset
        utils.visualization.create_2d_fields_from_pv_dataset(
            pv_dataset=pv_dataset,
            field_dir=field_dir,
            show_edges=show_edges,
            colormap=colormap,
            off_screen=off_screen,
            background_color=background_color,
            scalar_bar=scalar_bar,
            transparent_background=transparent_background,
            virtual_display=virtual_display)
