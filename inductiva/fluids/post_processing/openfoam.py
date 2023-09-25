"""Post-processing tools of fluid dynamics steady-state simulations.

This class implements various post-processing capabilities for
the visualizations associated. Namely:
    - Pressure over object;
    - Cutting plane;
    - Stream lines.

Current Support - OpenFOAM
"""
import os
from dataclasses import dataclass
from enum import Enum
from typing import Literal
import pathlib
import glob

from absl import logging

try:
    import pyvista as pv
except ImportError:
    pv = None

from inductiva import types, utils
from inductiva.utils import optional_deps


class SteadyStateOutput:
    """Post-Process steady-state simulation outputs.

    This class contains several methods to post-process the output 
    and visualize the results of a steady-state simulation where
    time-independent results are obtained.

    To be general we assume that a simulation is performed over a
    regular box - the domain - and a certain object placed inside.
    The object is either a small object inside the domain or a more
    general object that spreads all through the domain.
    """

    def __init__(self, sim_output_path: types.Path):
        """Initializes a `SteadyStateOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
            default_output_files_list: List of default output files
                that are computed on the backend and saved to the
                simulation output directory.
        """

        self.sim_output_path = sim_output_path
        self.full_output = self.inspect_output()

    def inspect_output(self):
        """Inspect the output of the simulation.
        
        Return True if the simulation output contains more files than
        the default ones. This is used to determine if the post-processing
        was computed on the backend or not."""

        # TODO: Transport these files in a simpler way from the task.
        # For such, a change to all of the post-processing output classes may
        # be required.
        default_output_files_list = [
            "pressure_field.vtk", "streamlines.vtk", "stdout.txt", "stderr.txt",
            "force_coefficients.csv", "xy_flow_slice.vtk", "xz_flow_slice.vtk",
            "yz_flow_slice.vtk", "object.obj"
        ]

        sim_output_files = glob.glob(os.path.join(self.sim_output_path, "**",
                                                  "*.*"),
                                     recursive=True)
        sim_output_files = [os.path.basename(file) for file in sim_output_files]

        return sorted(sim_output_files) != sorted(default_output_files_list)

    def get_last_iteration(self):
        """Get the last iteration of the simulation from simulation outputs."""

        # Sort all output folders and files by alphabetical order.
        outputs_dir_list = sorted(os.listdir(self.sim_output_path))
        # Second folder contains the last iteration outputs.
        last_iteration = float(outputs_dir_list[1])

        return last_iteration

    @optional_deps.needs_fluids_extra_deps
    def get_output_mesh(self):
        """Get domain and object mesh info at the steady-state."""

        if not self.full_output:
            object_mesh = pv.read(
                os.path.join(self.sim_output_path, "constant", "triSurface",
                             "object.obj"))
            return object_mesh

        # The OpenFOAM data reader from PyVista requires that a file named
        # "foam.foam" exists in the simulation output directory.
        # Create this file if it does not exist.
        foam_file_path = os.path.join(self.sim_output_path, "foam.foam")
        pathlib.Path(foam_file_path).touch(exist_ok=True)

        reader = pv.OpenFOAMReader(foam_file_path)
        reader.set_active_time_value(self.get_last_iteration())

        full_mesh = reader.read()
        domain_mesh = full_mesh["internalMesh"]
        object_mesh = full_mesh["boundary"]["object"]

        return domain_mesh, object_mesh

    @optional_deps.needs_fluids_extra_deps
    def get_object_pressure_field(self, save_path: types.Path = None):
        """Get a pressure scalar field over mesh points of the object.

        Returns:
            A MeshData object that allow to manipulate the data over a mesh
            and to render it.
        """

        if not self.full_output:
            logging.info("Pressure field was computed on the backend. "
                         "Args passed here were ignored.")
            object_mesh = pv.read(
                os.path.join(self.sim_output_path, "pressure_field.vtk"))
        else:
            logging.info("Fetching pressure field over the object.")
            _, object_mesh = self.get_output_mesh()

        field_notation = OpenFOAMPhysicalField["PRESSURE"].value
        pressure_field = MeshData(object_mesh, field_notation)

        if save_path is not None:
            save_path = utils.files.resolve_path(save_path)
            pressure_field.mesh.save(save_path)

        return pressure_field

    @optional_deps.needs_fluids_extra_deps
    def get_streamlines(self,
                        max_time: float = 100,
                        n_points: int = 100,
                        initial_step_length: float = 1,
                        source_radius: float = 0.7,
                        source_center: tuple = (0, 0, 0),
                        save_path: types.Path = None):
        """Get streamlines through the fluid/domain.
        
        The streamlines are obtained by seeding a set of points
        at the inlet of the WindTunnel.

        Args:
            simulation_time: Time value to obtain simulation mesh.
            max_time: Time used for integration of the streamlines.
                Not related with simulation time.
            n_points: Number of points to seed.
            initial_step_length: Initial step length for the streamlines.
            source_radius: Radius of the source of the streamlines.
            source_center: Center of the source of the streamlines.
            save_path: Path to save the streamlines. 
                Types of files permitted: .vtk, .ply, .stl
        """

        if not self.full_output:
            logging.info("Streamlines were computed on the backend. "
                         "Args passed here were ignored.")
            object_mesh = self.get_output_mesh()
            streamlines_mesh = pv.read(
                os.path.join(self.sim_output_path, "streamlines.vtk"))
        else:
            logging.info("Computing the streamlines on the fly.")
            mesh, object_mesh = self.get_output_mesh()

            streamlines_mesh = mesh.streamlines(
                max_time=max_time,
                n_points=n_points,
                initial_step_length=initial_step_length,
                source_radius=source_radius,
                source_center=source_center)

            if save_path is not None:
                save_path = utils.files.resolve_path(save_path)
                streamlines_mesh.save(save_path)

        return Streamlines(object_mesh, streamlines_mesh)

    @optional_deps.needs_fluids_extra_deps
    def get_flow_slice(self,
                       plane: Literal["xy", "xz", "yz"] = "xz",
                       origin: tuple = (0, 0, 0),
                       save_path: types.Path = None):
        """Get flow properties in a slice of the domain.
        
        Args:
            simulation_time: Time value to obtain simulation mesh.
            plane: Orientation of the plane to slice the domain.
            origin: Origin of the plane.
            save_path: Path to save the flow slice. 
                Types of files permitted: .vtk, .ply, .stl
        """

        if not self.full_output:
            logging.info("Flow slices were computed on the backend. "
                         "The origin and save_path arg are ignored.")
            object_mesh = self.get_output_mesh()
            flow_slice_name = plane + "_flow_slice.vtk"
            flow_slice = pv.read(
                os.path.join(self.sim_output_path, flow_slice_name))
        else:
            logging.info("Computing the flow on a slice.")
            mesh, object_mesh = self.get_output_mesh()

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
                save_path = utils.files.resolve_path(save_path)
                flow_slice.save(save_path)

        return FlowSlice(object_mesh, flow_slice)


@dataclass
class OpenFOAMPhysicalField(Enum):
    """Defines the notation used for physical field in OpenFOAM."""
    PRESSURE = "p"
    VELOCITY = "U"


class FlowSlice:
    """Render flow field in a plane of the domain."""

    def __init__(self, object_mesh, flow_slice):
        self.object_mesh = object_mesh
        self.mesh = flow_slice

    def render(self,
               physical_field: Literal["pressure", "velocity"] = "pressure",
               off_screen: bool = False,
               virtual_display: bool = False,
               background_color: str = "white",
               flow_cmap: str = "viridis",
               object_color: str = "white",
               save_path: types.Path = "flow_slice.png"):
        """Render flow property over domain."""

        if save_path is not None:
            save_path = utils.files.resolve_path(save_path)

        if virtual_display:
            # Set off_screen to True, since screens aren't available.
            off_screen = True
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

        plotter.add_mesh(self.object_mesh, color=object_color)
        plotter.add_mesh(self.mesh, scalars=field_notation, cmap=flow_cmap)
        plotter.reset_camera(bounds=self.mesh.bounds)
        plotter.show()
        if save_path is not None:
            plotter.screenshot(save_path, return_img=False)
            logging.info("Flow slice rendering saved to %s", save_path)
        plotter.close()


class Streamlines:
    """Class to render streamlines through the domain with object."""

    def __init__(self, object_mesh, streamlines):
        self.object_mesh = object_mesh
        self.mesh = streamlines

    def render(self,
               physical_field: Literal["pressure", "velocity"] = "pressure",
               off_screen: bool = False,
               virtual_display: bool = False,
               background_color: str = "white",
               flow_cmap: str = "viridis",
               view: Literal["isometric", "front", "rear", "top",
                             "side"] = "isometric",
               object_color: str = "white",
               streamline_radius: float = 0.1,
               save_path: types.Path = "streamlines.png"):
        """Render streamlines through domain."""

        if save_path is not None:
            save_path = utils.files.resolve_path(save_path)

        if virtual_display:
            # Set off_screen to True, since screens aren't available.
            off_screen = True
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

        plotter.add_mesh(self.mesh.tube(radius=streamline_radius),
                         scalars=field_notation,
                         cmap=flow_cmap)

        plotter.add_mesh(self.object_mesh, color=object_color)
        # Slide along the vectord defined from camera position to focal point,
        # until all of the meshes are visible.
        plotter.reset_camera(bounds=self.mesh.bounds, render=False)
        plotter.show()
        if save_path is not None:
            plotter.screenshot(save_path, return_img=False)
            logging.info("Streamlines rendering saved to %s", save_path)
        plotter.close()


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
            mesh: Mesh over the object obtained from the simulation.
        """
        self.mesh = pv.PolyData(mesh_data.points, faces=mesh_data.faces)
        self.scalar_name = scalar_name
        self.mesh.point_data[scalar_name] = mesh_data.point_data[scalar_name]
        self.mesh.cell_data[scalar_name] = mesh_data.cell_data[scalar_name]

    def render(self,
               off_screen: bool = False,
               background_color: str = "white",
               scalars_cmap: str = "viridis",
               virtual_display: bool = False,
               save_path: types.Path = "pressure_field.png"):
        """Render scalar field data over the mesh."""
        if save_path is not None:
            save_path = utils.files.resolve_path(save_path)

        if virtual_display:
            # Set off_screen to True, since screens aren't available.
            off_screen = True
            pv.start_xvfb()

        plotter = pv.Plotter(off_screen=off_screen)
        pv.global_theme.background = background_color
        plotter.add_mesh(self.mesh, scalars=self.scalar_name, cmap=scalars_cmap)
        plotter.show()
        if save_path is not None:
            plotter.screenshot(save_path, return_img=False)
            logging.info("Data rendering was saved to %s", save_path)
        plotter.close()
