"""Visualization processing of WindTunnel scenario.

This class implements various visualization capabilities for
the WindTunnel scenario. Namely:
    - Pressure over object;
    - Cutting plane;
    - StreamLines.

Currently, we only support the OpenFOAM simulator.
"""
import os
from dataclasses import dataclass
from enum import Enum
from typing import Literal
import csv

from base64 import b64encode
from IPython.display import HTML

import pyvista as pv

from inductiva.types import Path
from inductiva.utils.visualization import MeshData
from inductiva.utils import files


class WindTunnelSimulationOutput:
    """Post-Process WindTunnel simulation outputs.

    This class contains several methods to post-process the output 
    and visualize the results of a WindTunnel simulation.

    Current Support:
        OpenFOAM
    """

    def __init__(self, sim_output_path: Path, time_step: int = 100):
        """Initializes a `WindTunnelSimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
            time_step: Time step where we read the data (s).
                This time step needs to be a multiple of the
                output time step of the simulation.
        """

        self.sim_output_path = sim_output_path
        self.time_step = time_step
        self.object_data = self.get_object_data()

    def get_object_data(self):
        """Get aerodynamics measumerents at the object inside the WindTunnel.

        Current Support - OpenFOAM
        """

        reading_file = os.path.join(self.sim_output_path, "foam.foam")

        # Create reading file
        with open(reading_file, "w", encoding="utf-8") as file:
            file.close()

        # Initialize reader and define reading time-step
        reader = pv.POpenFOAMReader(reading_file)
        reader.set_active_time_value(self.time_step)

        mesh = reader.read()
        object_data = mesh["boundary"]["object"]

        return object_data

    def get_physical_field(self, physical_property: str = "pressure"):
        """Get a physical scalar field over mesh points for a certain time_step.

        Args:
            physical_property: Physical property to be read.
        Returns:
            A MeshData object that allow to manipulate the data over a mesh
            and to render it.
        """
        property_notation = OpenFOAMPhysicalProperty[
            physical_property.upper()].value
        physical_field = MeshData(self.object_data, property_notation)

        return physical_field

    def get_streamlines(self,
                        physical_property: Literal["pressure",
                                                   "velocity"] = "pressure"):
        """Get streamlines over the object in the WindTunnel.
        
        Args:
            physical_property: Physical property to be read.
        Returns:
            A MeshData object with the streamlines that represent the flow.
        """

        streamlines_path = os.path.join(self.sim_output_path, "postProcessing",
                                        "sets", "streamLines",
                                        str(self.time_step))

        property_notation = OpenFOAMPhysicalProperty[
            physical_property.upper()].value
        streamlines_file = "track0_" + property_notation + ".vtk"

        streamlines_mesh = pv.read(
            os.path.join(streamlines_path, streamlines_file))

        return streamlines_mesh

    def get_flow_plane(self,
                       physical_property: Literal["pressure",
                                                  "velocity"] = "pressure"):
        """Get flow properties in a plane of the domain in WindTunnel.
        
        Args:
            physical_property: Physical property to be read.
        Returns:
            A MeshData object with aerodynamic measurements over a plane."""

        cutting_plane_path = os.path.join(self.sim_output_path,
                                          "postProcessing", "cuttingPlane",
                                          str(self.time_step))

        # Obtain notation for the physical property for the simulator.
        property_notation = OpenFOAMPhysicalProperty[
            physical_property.upper()].value
        cutting_plane_file = property_notation + "_yNormal.vtk"

        cutting_plane_mesh = pv.read(
            os.path.join(cutting_plane_path, cutting_plane_file))

        return cutting_plane_mesh

    def get_force_coefficients(self, save_path: Path = None):
        """Get the force coefficients of the object in the WindTunnel.
        
        The force coefficients are provided in a .dat file during the
        simulation run-time. This file contains 8 lines that are provide
        the general input information. In this function, we read the file,
        ignore the first 8 lines and read the force coefficients for the 
        time_step chosen.

        Args:
            save_path: Path to save the force coefficients in a .csv file.
        """

        num_header_lines = 8
        force_coefficients_path = os.path.join(self.sim_output_path,
                                               "postProcessing", "forceCoeffs1",
                                               "0", "forceCoeffs.dat")
        force_coefficients = []

        with open(force_coefficients_path, "r",
                  encoding="utf-8") as forces_file:
            for index, line in enumerate(forces_file.readlines()):
                # Pick the line 8 of the file:
                # [#, Time, Cm, Cd, Cl, Cl(f), Cl(r)] and remove the # column
                if index == num_header_lines:
                    force_coefficients.append(line.split()[1:])
                # Add the force coefficients for the time_step chosen
                elif index == num_header_lines + self.time_step + 1:
                    force_coefficients.append(line.split())

        if save_path:
            with open(save_path, "w", encoding="utf-8") as csv_file:
                csv_writer = csv.writer(csv_file)
                csv_writer.writerows(force_coefficients)

        return force_coefficients

    def render_flow(self,
                    flow_property_mesh,
                    physical_property: Literal["pressure",
                                               "velocity"] = "pressure",
                    virtual_display: bool = True,
                    background_color: str = "black",
                    flow_cmap: str = "viridis",
                    object_color: str = "white",
                    save_path: Path = None):
        """Render flow property over the object in the WindTunnel.
        
        Args:
            flow_property_mesh: MeshData object with the flow property.
            physical_property: Physical property to be read.
            virtual_display: Uses a virtual display to render the plot.
                Essential, to render in a remote server.
            background_color: Background color of the plot.
            flow_cmap: Colormap for the flow property.
            object_color: Color of the object.
            save_path: Path to save the plot. If None, the
                plot is not saved.
        """
        if save_path is not None:
            save_path = files.resolve_path(save_path)

        off_screen = False

        if virtual_display:
            off_screen = True
            pv.start_xvfb()

        # Obtain notation for the physical property for the simulator.
        property_notation = OpenFOAMPhysicalProperty[
            physical_property.upper()].value

        plotter = pv.Plotter(off_screen=off_screen)
        plotter.background_color = background_color
        plotter.add_mesh(self.object_data, color=object_color)
        plotter.add_mesh(flow_property_mesh,
                         scalars=property_notation,
                         cmap=flow_cmap)
        plotter.view_xz()
        plotter.show(screenshot=save_path)
        plotter.close()

        with open(save_path, "rb") as file_path:
            png = file_path.read()
        png_url = "data:image/png;base64," + b64encode(png).decode()

        return HTML(f"""
                <img src="{png_url}" type="image/png" width="600">
        """)


@dataclass
class OpenFOAMPhysicalProperty(Enum):
    """Notation for physical properties in OpenFOAM."""
    PRESSURE = "p"
    VELOCITY = "U"
