"""Post-processing of fluid dynamics steady-state simulations.

This class implements various post-processing capabilities with
visuals associated. Namely:
    - Pressure over object;
    - Cutting plane;
    - StreamLines.
"""
import os
from dataclasses import dataclass
import enum
import typing
import pathlib

import pyvista as pv

import inductiva


class SteadyStateOutput:
    """Post-Process steady-state simulation outputs.

    This class contains several methods to post-process the output 
    and visualize the results of steady-state simulations where
    static simulations are performed.

    To be general we assume that the simulation is performed in a
    box and a certain object is placed inside the box. The object
    can either refer to a small object inside the domain or a more
    general object that spreads all through the domain.
    """

    def __init__(self, sim_output_path: inductiva.types.Path):
        """Initializes a `SteadyStateOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
        """

        self.sim_output_path = sim_output_path

    def get_mesh_at_time(self, simulation_time: float):  # pylint: disable=unused-argument
        """Get domain and object mesh info after steady-state simulation.

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

    def get_object_pressure_field(self,
                                  simulation_time: float = 50,
                                  save_path: inductiva.types.Path = None):
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

        field_notation = OpenFOAMPhysicalField["PRESSURE"].value
        physical_field = MeshData(object_mesh, field_notation)

        if save_path is not None:
            save_path = inductiva.utils.files.resolve_path(save_path)
            physical_field.mesh.save(save_path)

        return physical_field

    def get_streamlines(self,
                        simulation_time: float = 50,
                        max_time: float = 100,
                        n_points: int = 100,
                        initial_step_length: float = 1,
                        source_radius: float = 0.7,
                        save_path: inductiva.types.Path = None):
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

        domain_mesh, object_mesh = self.get_mesh_at_time(simulation_time)

        inlet_position = (domain_mesh.bounds[0], 0, 1)

        streamlines_mesh = domain_mesh.streamlines(
            max_time=max_time,
            n_points=n_points,
            initial_step_length=initial_step_length,
            source_radius=source_radius,
            source_center=inlet_position)

        if save_path is not None:
            save_path = inductiva.utils.files.resolve_path(save_path)
            streamlines_mesh.save(save_path)

        return Streamlines(streamlines_mesh, object_mesh)

    def get_flow_slice(self,
                       simulation_time: float = 50,
                       plane: typing.Literal["xy", "xz", "yz"] = "xz",
                       origin: tuple = (0, 0, 0),
                       save_path: inductiva.types.Path = None):
        """Get flow properties in a slice of the domain in WindTunnel.
        
        Args:
            simulation_time: Time value to obtain simulation mesh.
            plane: Orientation of the plane to slice the domain.
            origin: Origin of the plane.
            save_path: Path to save the flow slice. 
                Types of files permitted: .vtk, .ply, .stl
        """

        domain_mesh, object_mesh = self.get_mesh_at_time(simulation_time)

        if plane == "xy":
            normal = (0, 0, 1)
        elif plane == "yz":
            normal = (1, 0, 0)
        elif plane == "xz":
            normal = (0, 1, 0)
        else:
            raise ValueError("Invalid view.")

        flow_slice = domain_mesh.slice(normal=normal, origin=origin)

        if save_path is not None:
            save_path = inductiva.utils.files.resolve_path(save_path)
            flow_slice.save(save_path)

        return FlowSlice(flow_slice, object_mesh)


@dataclass
class OpenFOAMPhysicalField(enum.Enum):
    """Defines the notation used for physical field in OpenFOAM."""
    PRESSURE = "p"
    VELOCITY = "U"


class MeshData:
    """MeshData class that allows for mesh data manipulation and render.

    Process outputs of simulation over meshes. We assume that scalar
    fields are defined over the points of the mesh. In this way, we have
    a general method to manipulate the data and render the specification
    data. As of now, this serves a niche, but the goal
    is to integrate this more generally later on.

    Example:
    Imagine we run a WindTunnel simulation with an object inside
    for which a mesh was generated. A possible output of the simulation is
    a mesh of the object with scalar fields over the points defining a certain
    physical property, e.g., pressure.

    Assume this output mesh is `object_mesh`. To process the pressure field over
    the data, we can do `pressure_field = MeshData(object_mesh, "p")`

    `pressure_field.mesh` contains a mesh structure with the pressure field over
    the mesh.

    `pressure_field.render()` renders the pressure field property over the mesh.
    """

    def __init__(self, mesh_data, scalar_name: str = None):
        """Initialize a `MeshData` type object.

        Args:
            mesh_data: pyvista mesh (PolyData or Unstructured) which
                contains all of the simulated outputs over the mesh.
            scalar_name: string that defines the scalar we want to
                manipulate and render. This names depends on the mesh_data
                array names and on the specific simulator used.

        Attributes:
            mesh: mesh over the object obtained from
        """
        self.mesh = pv.PolyData(mesh_data.points, faces=mesh_data.faces)
        self.scalar_name = scalar_name
        self.mesh.point_data[scalar_name] = mesh_data.point_data[scalar_name]
        self.mesh.cell_data[scalar_name] = mesh_data.cell_data[scalar_name]

    def render(self,
               off_screen: bool = False,
               background_color: str = "black",
               colormap: str = "viridis",
               virtual_display: bool = False,
               save_path: inductiva.types.Path = None):
        """Render scalar field data over the mesh.
        
        Args:
            off_screen: Renders the plot off the screen and avoids the
                deployment of a visualization widget. To obtain the results
                from this, one needs to give a path to save the frame.
            background_color: Color of the background of the plot.
            colormap: Colormap to use for the scalar field being plotted.
            virtual_display: If True, a virtual display is created. Essential
            to render on a server.
            save_path: Path to save the frame. Types of files permitted: .png
        """
        if save_path is not None:
            save_path = inductiva.utils.files.resolve_path(save_path)

        if virtual_display:
            pv.start_xvfb()

        plotter = pv.Plotter(off_screen=off_screen)
        pv.global_theme.background = background_color
        plotter.add_mesh(self.mesh, scalars=self.scalar_name, cmap=colormap)
        plotter.show(screenshot=save_path)
        plotter.close()


class FlowSlice:
    """Render flow field in a plane of the domain in WindTunnel."""

    def __init__(self, flow_slice, object_mesh):
        self.flow_mesh = flow_slice
        self.object_mesh = object_mesh

    def render_frame(self,
                     physical_field: typing.Literal["pressure",
                                                    "velocity"] = "pressure",
                     off_screen: bool = False,
                     virtual_display: bool = False,
                     background_color: str = "black",
                     colormap: str = "viridis",
                     object_color: str = "white",
                     save_path: inductiva.types.Path = None):
        """Render flow property over the object in the WindTunnel.
        
        Args:
            physical_field: Physical field to render. Options are
                pressure and velocity.
            off_screen: Renders the plot off the screen and avoids the
                deployment of a visualization widget. To obtain the results
                from this, one needs to give a path to save the frame.
            background_color: Color of the background of the plot.
            colormap: Colormap to use for the scalar field.
            virtual_display: If True, a virtual display is created. Essential
            to render on a server.
            object_color: Color of the object mesh.
            save_path: Path to save the frame. Types of files permitted: .png
        """

        if save_path is not None:
            save_path = inductiva.utils.files.resolve_path(save_path)

        if virtual_display:
            pv.start_xvfb()

        plotter = pv.Plotter(off_screen=off_screen)
        plotter.background_color = background_color

        center = self.flow_mesh.center
        normal = -self.flow_mesh.compute_normals()["Normals"].mean(axis=0)
        plotter.camera.position = center + normal
        plotter.camera.focal_point = center
        plotter.camera.zoom(1.2)

        # Obtain notation for the physical field for the simulator.
        field_notation = OpenFOAMPhysicalField[physical_field.upper()].value

        plotter.add_mesh(self.object_mesh, color=object_color)
        plotter.add_mesh(self.flow_mesh, scalars=field_notation, cmap=colormap)
        plotter.reset_camera(bounds=self.flow_mesh.bounds)
        plotter.show(screenshot=save_path)
        plotter.close()


class Streamlines:
    """Class to render streamlines over the object in the WindTunnel."""

    def __init__(self, streamlines, object_mesh):
        self.streamlines_mesh = streamlines
        self.object_mesh = object_mesh

    def render_frame(self,
                     physical_field: typing.Literal["pressure",
                                                    "velocity"] = "pressure",
                     off_screen: bool = False,
                     virtual_display: bool = False,
                     background_color: str = "black",
                     colormap: str = "viridis",
                     view: typing.Literal["isometric", "front", "rear", "top",
                                          "side"] = "isometric",
                     object_color: str = "white",
                     save_path: inductiva.types.Path = None):
        """Render streamlines over the object in the WindTunnel.
        
        Args:
            physical_field: Physical field to render. Options are
                pressure and velocity.
            off_screen: Renders the plot off the screen and avoids the
                deployment of a visualization widget. To obtain the results
                from this, one needs to give a path to save the frame.
            background_color: Color of the background of the plot.
            colormap: Colormap to use for the scalar field for the streamlines
                mesh.
            virtual_display: If True, a virtual display is created. Essential
            to render on a server.
            view: Camera view of the plot. Options are isometric, front, rear,
                top, side.
            object_color: Color of the object mesh.
            save_path: Path to save the frame. Types of files permitted: .png
        """

        if save_path is not None:
            save_path = inductiva.utils.files.resolve_path(save_path)

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

        plotter.add_mesh(self.streamlines_mesh.tube(radius=0.01),
                         scalars=field_notation,
                         cmap=colormap)
        plotter.add_mesh(self.object_mesh, color=object_color)
        # Slide along the vector defined from camera position to focal point,
        # until all of the meshes are visible.
        plotter.reset_camera(bounds=self.streamlines_mesh.bounds, render=False)
        plotter.show(screenshot=save_path)
        plotter.close()
