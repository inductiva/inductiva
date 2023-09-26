"""Utils for structures case."""


def remove_time_values_from_xdmf(input_file_path: str,
                                 output_file_path: str) -> None:
    """Removes Time Value lines from an XDMF file.

    This method is crucial due to PyVista's limitations in inspecting time
    values directly. It reads an XDMF file, removes lines containing Time 
    Values, and writes the modified content to a new file.

    Args:
        input_file_path (str): The path to the input XDMF file.
        output_file_path (str): The path to save the modified XDMF content.
    """

    # Read the XDMF file
    with open(input_file_path, "r", encoding="utf-8") as file:
        xdmf_content = file.readlines()

    # Create an empty list to store the modified content
    new_xdmf_content = []

    # Iterate through the lines in the file
    for line in xdmf_content:

        # Check if the line contains a Time Value
        if "<Time Value=" not in line:

            # If the line does not contain a Time Value, add it to the new file
            new_xdmf_content.append(line)

    # Write the modified content back to the output file
    with open(output_file_path, "w", encoding="utf-8") as file:
        file.writelines(new_xdmf_content)
