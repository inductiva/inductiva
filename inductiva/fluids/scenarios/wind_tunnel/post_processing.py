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
import pathlib

import pyvista as pv

from inductiva.types import Path
from inductiva.utils.visualization import MeshData
from inductiva.utils import files


class WindTunnelOutput:
    """Post-Process WindTunnel simulation outputs.

    Current Support:
        OpenFOAM
    """

    def __init__(self, sim_output_path: Path):
        """Initializes a `WindTunnelSimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
        """

        self.sim_output_path = sim_output_path

    def get_mesh_data(self, time_step):  # pylint: disable=unused-argument
        """Get aerodynamics data over an object inside the WindTunnel.

        Current Support - OpenFOAM
        """

        # The OpenFOAM data reader from PyVista requires that a file named
        # "foam.foam" exists in the simulation output directory.
        # Create this file if it does not exist.
        foam_file_path = os.path.join(self.sim_output_path, "foam.foam")
        pathlib.Path(foam_file_path).touch(exist_ok=True)

        reader = pv.OpenFOAMReader(foam_file_path)
        reading_file = os.path.join(self.sim_output_path, "foam.foam")

        # Initialize reader and define reading time-step
        reader = pv.POpenFOAMReader(reading_file)
        reader.set_active_time_value(time_step)

        full_mesh = reader.read()
        domain_mesh = full_mesh["internalMesh"]
        object_mesh = full_mesh["boundary"]["object"]

        return domain_mesh, object_mesh

    def get_physical_field(self,
                           physical_property: str = "pressure",
                           time_step: float = 50):
        """Get a physical scalar field over mesh points for a certain time_step.

        Returns:
            A MeshData object that allow to manipulate the data over a mesh
            and to render it.
        """

        _, object_mesh = self.get_mesh_data(time_step)

        property_notation = OpenFOAMPhysicalProperty[
            physical_property.upper()].value
        physical_field = MeshData(object_mesh, property_notation)

        return physical_field

    def get_streamlines(self,
                        time_step: float = 50,
                        max_time: float = 100,
                        n_points: int = 100,
                        initial_step_length: float = 1,
                        source_radius: float = 0.7,
                        save_path: Path = None):
        """Get streamlines over the object in the WindTunnel.
        
        The streamlines are obtained by seeding a set of points
        at the inlet of the WindTunnel.

        Args:
            time_step: Simulation time step to obtain the streamlines.
            max_time: Maximum time for the streamlines.
            n_points: Number of points to seed.
            initial_step_length: Initial step length for the streamlines.
            source_radius: Radius of the source of the streamlines.
        """

        mesh, _ = self.get_mesh_data(time_step)

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
                       time_step: float = 50,
                       plane: Literal["xy", "xz", "yz"] = "xy",
                       origin: tuple = (0, 0, 0),
                       save_path: Path = None):
        """Get flow properties in a plane of the domain in WindTunnel."""

        mesh, _ = self.get_mesh_data(time_step)

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


@dataclass
class OpenFOAMPhysicalProperty(Enum):
    """Defines the notation used for physical properties in OpenFOAM."""
    PRESSURE = "p"
    VELOCITY = "U"


class FlowSlice:
    """Render flow properties in a plane of the domain in WindTunnel."""

    def __init__(self, flow_slice):
        self.mesh = flow_slice

    def render_frame(self,
                     object_mesh: pv.PolyData = None,
                     physical_property: Literal["pressure",
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
        normal = self.mesh.compute_normals()["Normals"].mean(axis=0)
        plotter.camera.position = center + normal
        plotter.camera.focal_point = center
        plotter.camera.zoom(1.2)

        # Obtain notation for the physical property for the simulator.
        property_notation = OpenFOAMPhysicalProperty[
            physical_property.upper()].value

        plotter.add_mesh(object_mesh, color=object_color)
        plotter.add_mesh(self.mesh, scalars=property_notation, cmap=flow_cmap)
        plotter.reset_camera(bounds=self.mesh.bounds)
        plotter.show(screenshot=save_path)
        plotter.close()


class Streamlines:
    """Class to render streamlines over the object in the WindTunnel."""

    def __init__(self, streamlines):
        self.mesh = streamlines

    def render_frame(self,
                     object_mesh: pv.PolyData = None,
                     physical_property: Literal["pressure",
                                                "velocity"] = "pressure",
                     off_screen: bool = False,
                     virtual_display: bool = True,
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

        # Obtain notation for the physical property for the simulator.
        property_notation = OpenFOAMPhysicalProperty[
            physical_property.upper()].value

        plotter.add_mesh(self.mesh.tube(radius=0.01),
                         scalars=property_notation,
                         cmap=flow_cmap)
        plotter.add_mesh(object_mesh, color=object_color)
        # Slide along the vectord defined from camera position to focal point,
        # until all of the meshes are visible.
        plotter.reset_camera(bounds=self.mesh.bounds, render=False)
        plotter.show(screenshot=save_path)
        plotter.close()
