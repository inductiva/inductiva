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
import pathlib

import pyvista as pv

from inductiva.types import Path
from inductiva.utils.visualization import MeshData
from inductiva.utils import files


class WindTunnelOutput:
    """Post-Process WindTunnel simulation outputs.

    This class contains several methods to post-process the output 
    and visualize the results of a WindTunnel simulation.

    Current Support:
        OpenFOAM
    """

    def __init__(self, sim_output_path: Path):
        """Initializes a `WindTunnelSimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
        """

        self.sim_output_path = sim_output_path

    def get_mesh_at_time(self, simulation_time):  # pylint: disable=unused-argument
        """Get domain and object mesh info after WindTunnel simulation.

        Current Support - OpenFOAM

        Args:
            simulation_time: Time value to obtain simulation mesh.

        """

        # The OpenFOAM data reader from PyVista requires that a file named
        # "foam.foam" exists in the simulation output directory.
        # Create this file if it does not exist.
        foam_file_path = os.path.join(self.sim_output_path, "foam.foam")
        pathlib.Path(foam_file_path).touch(exist_ok=True)

        reader = pv.OpenFOAMReader(foam_file_path)
        reader.set_active_time_value(simulation_time)

        full_mesh = reader.read()
        domain_mesh = full_mesh["internalMesh"]
        object_mesh = full_mesh["boundary"]["object"]

        return domain_mesh, object_mesh

    def get_object_physical_field(self,
                                  physical_field: str = "pressure",
                                  simulation_time: float = 50,
                                  save_path: Path = None):
        """Get a physical scalar field over mesh points.

        Args:
            simulation_time: Time value to obtain simulation mesh.

        Args:
            physical_property: Physical property to be read.
        Returns:
            A MeshData object that allow to manipulate the data over a mesh
            and to render it.
        """

        _, object_mesh = self.get_mesh_at_time(simulation_time)

        field_notation = OpenFOAMPhysicalField[physical_field.upper()].value
        physical_field = MeshData(object_mesh, field_notation)

        if save_path is not None:
            save_path = files.resolve_path(save_path)
            physical_field.mesh.save(save_path)

        return physical_field

    def get_streamlines(self,
                        simulation_time: float = 50,
                        max_time: float = 100,
                        n_points: int = 100,
                        initial_step_length: float = 1,
                        source_radius: float = 0.7,
                        save_path: Path = None):
        """Get streamlines through the fluid/domain in the WindTunnel.
        
        The streamlines are obtained by seeding a set of points
        at the inlet of the WindTunnel.

        Args:
            simulation_time: Time value to obtain simulation mesh.
            max_time: Time used for integration of the streamlines.
                Not related with simulation time.
            n_points: Number of points to seed.
            initial_step_length: Initial step length for the streamlines.
            source_radius: Radius of the source of the streamlines.
            save_path: Path to save the streamlines. 
                Types of files permitted: .vtk, .ply, .stl
        """

        mesh, _ = self.get_mesh_at_time(simulation_time)

        inlet_position = (mesh.bounds[0], 0, 1)

        streamlines_mesh = mesh.streamlines(
            max_time=max_time,
            n_points=n_points,
            initial_step_length=initial_step_length,
            source_radius=source_radius,
            source_center=inlet_position)

        if save_path is not None:
            save_path = files.resolve_path(save_path)
            streamlines_mesh.save(save_path)

        return Streamlines(streamlines_mesh)

    def get_flow_slice(self,
                       simulation_time: float = 50,
                       plane: Literal["xy", "xz", "yz"] = "xz",
                       origin: tuple = (0, 0, 0),
                       save_path: Path = None):
        """Get flow properties in a slice of the domain in WindTunnel.
        
        Args:
            simulation_time: Time value to obtain simulation mesh.
            plane: Orientation of the plane to slice the domain.
            origin: Origin of the plane.
            save_path: Path to save the flow slice. 
                Types of files permitted: .vtk, .ply, .stl
        """

        mesh, _ = self.get_mesh_at_time(simulation_time)

        if plane == "xy":
            normal = (0, 0, 1)
        elif plane == "yz":
            normal = (1, 0, 0)
        elif plane == "xz":
            normal = (0, 1, 0)
        else:
            raise ValueError("Invalid view.")

        flow_slice = mesh.slice(normal=normal, origin=origin)

        if save_path is not None:
            save_path = files.resolve_path(save_path)
            flow_slice.save(save_path)

        return FlowSlice(flow_slice)

    def get_force_coefficients(self,
                               simulation_time: float = 50,
                               save_path: Path = None):
        """Get the force coefficients of the object in the WindTunnel.
        
        The force coefficients are provided in a .dat file during the
        simulation run-time. This file contains 8 lines that are provide
        the general input information. In this function, we read the file,
        ignore the first 8 lines and read the force coefficients for the 
        simulation_time chosen.

        Args:
            simulation_time: Time value to obtain simulation mesh.
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
                # Add the force coefficients for the simulation time chosen
                elif index == num_header_lines + simulation_time + 1:
                    force_coefficients.append(line.split())

        if save_path:
            with open(save_path, "w", encoding="utf-8") as csv_file:
                csv_writer = csv.writer(csv_file)
                csv_writer.writerows(force_coefficients)

        return force_coefficients


@dataclass
class OpenFOAMPhysicalField(Enum):
    """Defines the notation used for physical field in OpenFOAM."""
    PRESSURE = "p"
    VELOCITY = "U"


class FlowSlice:
    """Render flow field in a plane of the domain in WindTunnel."""

    def __init__(self, flow_slice):
        self.mesh = flow_slice

    def render_frame(self,
                     object_mesh: pv.PolyData = None,
                     physical_field: Literal["pressure",
                                             "velocity"] = "pressure",
                     off_screen: bool = False,
                     virtual_display: bool = False,
                     background_color: str = "black",
                     flow_cmap: str = "viridis",
                     object_color: str = "white",
                     save_path: Path = None):
        """Render flow property over the object in the WindTunnel."""

        if save_path is not None:
            save_path = files.resolve_path(save_path)

        if virtual_display:
            pv.start_xvfb()

        plotter = pv.Plotter(off_screen=off_screen)
        plotter.background_color = background_color

        center = self.mesh.center
        normal = -self.mesh.compute_normals()["Normals"].mean(axis=0)
        plotter.camera.position = center + normal
        plotter.camera.focal_point = center
        plotter.camera.zoom(1.2)

        # Obtain notation for the physical field for the simulator.
        field_notation = OpenFOAMPhysicalField[physical_field.upper()].value

        if object_mesh:
            plotter.add_mesh(object_mesh, color=object_color)
        plotter.add_mesh(self.mesh, scalars=field_notation, cmap=flow_cmap)
        plotter.reset_camera(bounds=self.mesh.bounds)
        plotter.show(screenshot=save_path)
        plotter.close()


class Streamlines:
    """Class to render streamlines over the object in the WindTunnel."""

    def __init__(self, streamlines):
        self.mesh = streamlines

    def render_frame(self,
                     object_mesh: pv.PolyData = None,
                     physical_field: Literal["pressure",
                                             "velocity"] = "pressure",
                     off_screen: bool = False,
                     virtual_display: bool = False,
                     background_color: str = "black",
                     flow_cmap: str = "viridis",
                     view: Literal["isometric", "front", "rear", "top",
                                   "side"] = "isometric",
                     object_color: str = "white",
                     save_path: Path = None):
        """Render streamlines over the object in the WindTunnel."""

        if save_path is not None:
            save_path = files.resolve_path(save_path)

        if virtual_display:
            pv.start_xvfb()

        plotter = pv.Plotter(off_screen=off_screen)
        plotter.background_color = background_color

        # Set camera position for a nice view.
        if view == "isometric":
            plotter.view_vector([-1, -2, 1], viewup=[0.31, 0.90, 0.29])
        elif view == "front":
            plotter.view_yz(negative=True)
        elif view == "rear":
            plotter.view_yz()
        elif view == "top":
            plotter.view_xy()
        elif view == "side":
            plotter.view_xz()
        else:
            raise ValueError("Invalid view.")

        # Obtain notation for the physical field for the simulator.
        field_notation = OpenFOAMPhysicalField[physical_field.upper()].value

        plotter.add_mesh(self.mesh.tube(radius=0.01),
                         scalars=field_notation,
                         cmap=flow_cmap)
        if object_mesh:
            plotter.add_mesh(object_mesh, color=object_color)
        # Slide along the vectord defined from camera position to focal point,
        # until all of the meshes are visible.
        plotter.reset_camera(bounds=self.mesh.bounds, render=False)
        plotter.show(screenshot=save_path)
        plotter.close()
