"""Post-process FEniCSx simulation outputs."""

import os

import pyvista as pv

import inductiva


class FEniCSxSimulationOutput:
    """Post-process FEniCSx simulation outputs."""

    def __init__(self, sim_output_path: inductiva.types.Path):
        """Initializes a 'FEniCSxSimulationOutput' object.

        Args:
            sim_output_path: Path to simulation output files.
            """
        self.sim_output_path = sim_output_path

    def remove_time_value_lines(self, output_file_path):
        # Read the Xdmf file
        with open(self.sim_output_path, 'r') as file:
            xdmf_content = file.readlines()

        # Create an empty list to store the modified content
        new_xdmf_content = []

        # Iterate through the lines in the file
        for line in xdmf_content:
            # Check if the line contains a Time Value
            if not '<Time Value=' in line:
                new_xdmf_content.append(line)

        # Write the modified content back to the output file
        with open(output_file_path, 'w') as file:
            file.writelines(new_xdmf_content)

    @inductiva.utils.optional_deps.needs_structures_extra_deps
    def render(self,
               field_names: list = None,
               scalar_bar: bool = True,
               transparent_background: bool = True) -> None:
        """Render the simulation fields as images.

        Args:
            field_names (list, optional): List of field names to render. Default
              is None.
            scalar_bar (bool, optional): Whether to include a scalar bar legend
              in the visualization. Default is True.
            transparent_background (bool, optional): Whether the background of
              the saved images should be transparent. Default is True.
        """

        # Read simulation output data
        reader = pv.get_reader(self.sim_output_path)

        data = reader.read()

        if field_names is None:

            # Default field names to render
            field_names = [
                "displacement_x", "displacement_y", "displacement_z",
                "stress_xx", "stress_xy", "stress_xz", "stress_yx", "stress_yy",
                "stress_yz", "stress_zx", "stress_zy", "stress_zz", "strain_xx",
                "strain_xy", "strain_xz", "strain_yx", "strain_yy", "strain_yz",
                "strain_zx", "strain_zy", "strain_zz", "von_mises"
            ]

        # Get the list of available keys (field names) in the PyVista dataset
        keys = list(data.keys())

        # Iterate through the specified field names to render
        for field_name in field_names:
            if field_name != "von_mises":

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

                # Extract the specific subfield array from the dataset
                field_array = data[keys[field_id]].get_array(0)[:, subfield_id]

            elif field_name == "von_mises":
                # For the special case of "von_mises" field, find its index in
                # the dataset keys
                field_id = keys.index(field_name)

                # Extract the "von_mises" field array from the dataset
                field_array = data[keys[field_id]].get_array(0)

            # Define the path to save the field visualization as an image
            field_path = os.path.join(self.sim_output_path, "fields",
                                      f"{field_name}.png")

            # Create a 2D field visualization from the XDMF data
            inductiva.utils.visualization.create_2d_field_from_xdmf(
                xdmf_path=self.sim_output_path,
                field_path=field_path,
                field_array=field_array,
                file_name=field_name,
                scalar_bar=scalar_bar,
                transparent_background=transparent_background)
