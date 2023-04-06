"""Post-processing utilities for .vtk data"""

import os
import pathlib


def get_sorted_vtk_files(data_dir: str):
    """Returns a list of sorted vtk files in a directory.
    
    Order a set of .vtk files of the form
    ['name_1.vtk', 'name_2.vtk',...,'name_10.vtk',...,'name_n.vtk'].
    The sorting methods for list, list.sort() or sorted(list), set
    'name_10.vtk' before 'name_2.vtk', which is not representative
    of the time series.

    In this function we provide a correct method for this sorting
    """

    if not os.path.exists(data_dir):
        raise IOError(f"Directory '{data_dir}' does not exist.")

    # Get a list of the files in the data directory.
    files = os.scandir(data_dir)

    # The files must have .vtk extension.
    files = [file for file in files if pathlib.Path(file.path).suffix == ".vtk"]

    # Sort the files to be read according to [file_key].
    def get_alphanum_key(file):
        file_name = pathlib.Path(file.path).stem
        file_name_splits = file_name.split("_")
        file_key = file_name_splits[-1]
        return int(file_key)

    files = sorted(files, key=get_alphanum_key)
    return files
