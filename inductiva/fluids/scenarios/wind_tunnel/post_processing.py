"""Visualization processing of WindTunnel scenario.

This class implements various visualization capabilities for
the WindTunnel scenario. Namely:
    - Pressure over object; 
    - Cutting plane;
    - StreamLines.

Currently, we only support the OpenFOAM simulator.
"""
import os
import pyvista as pv

from inductiva.types import Path
from inductiva.utils.visualization import MeshData


class WindTunnelSimulationOutput:
    """Post process WindTunnel simulation outputs.
    
    Current Support:
        OpenFOAM
    """

    def __init__(self,
                 sim_output_path: Path,
                 time_step: int):
        """Initializes an `WindTunnelSimulationOutput` object.

        Args:
            simulator: Simulator object.
            sim_output_path: Path to simulation output files.
            time_step: Time step where we read the data.
        
        Attributes:
            sim_output_dir: path to the simulation directory
            object_path: path to the object insert in the WindTunnel
            post_processing_path: path to the post-processing data.
        """

        self.sim_output_dir = sim_output_path
        self.object_data = self.get_object_data(time_step)

    def get_object_data(self, time_step: int): # pylint: disable=unused-argument
        """Get aerodynamics data over an object set inside the WindTunnel
        
        Current Support - OpenFOAM

        Args:
            time_step: Time step for which we want to visualize
                aerodynamics data.
        """

        reading_file = os.path.join(self.sim_output_dir, "foam.foam")

        # Create reading file
        with open(reading_file, "w", encoding="utf-8") as file:
            file.close()

        # Initialize reader and define reading time-step
        reader = pv.POpenFOAMReader(reading_file)
        reader.set_active_time_value(time_step)

        mesh = reader.read()
        object_data = mesh["boundary"]["object"]

        return object_data

    def get_pressure_field(self): # pylint: disable=unused-argument
        """Get pressure field over mesh points for a certain time_step.
        
        Returns:
            A MeshData object that allow to manipulate data over a mesh
            and to render it.
        """
        pressure_field = MeshData(self.object_data, "p")

        return pressure_field

    def get_streamlines(self,
                        time_step: int,
                        physical_property: str = "pressure"):
        """Get streamlines over the object we insert in the WindTunnel."""

        streamlines_path = os.path.join(
            self.sim_output_dir, "postProcessing",
            "sets", "streamLines",
            str(time_step))

        if physical_property.lower() == "pressure":
            streamlines_mesh = pv.read(
            os.path.join(streamlines_path, "track0_p.vtk"))
        elif physical_property.lower() == "velocity":
            streamlines_mesh = pv.read(
            os.path.join(streamlines_path, "track0_U.vtk"))
        else:
            raise ValueError(f"Physical property `{physical_property}` not available "
                             "or supported")

        return streamlines_mesh

    def get_flow_plane(self,
                        time_step: int,
                        physical_property: str = "pressure"):
        """Get flow properties in a plane of the domain in WindTunnel."""
        cutting_plane_path = os.path.join(
            self.sim_output_dir, "postProcessing",
            "cuttingPlane", str(time_step))

        if physical_property.lower() == "pressure":
            cutting_plane_mesh = pv.read(os.path.join(
                cutting_plane_path, "p_yNormal.vtk"))
        elif physical_property.lower() == "velocity":
            cutting_plane_mesh = pv.read(os.path.join(
                cutting_plane_path, "U_yNormal.vtk"))
        else:
            raise ValueError(f"Physical property `{physical_property}` not available "
                             "or supported")

        return cutting_plane_mesh


    def render_flow(self, flow_property_mesh,
                    physical_property: str = "pressure"):
        """Render flow property over the object we insert in the WindTunnel."""

        plotter = pv.Plotter()
        plotter.add_mesh(self.object_data)
        if physical_property.lower() == "pressure":
            plotter.add_mesh(flow_property_mesh, scalars="p")
        elif physical_property.lower() == "velocity":
            plotter.add_mesh(flow_property_mesh, scalars="U")
        else:
            raise ValueError(f"Physical property `{physical_property}` not available "
                             "or supported")
        plotter.show()
        plotter.close()
