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
    """Post-Process WindTunnel simulation outputs.
    
    Current Support:
        OpenFOAM
    """

    def __init__(self, sim_output_path: Path, time_step: int):
        """Initializes an `WindTunnelSimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
            time_step: Time step where we read the data.
        """

        self.sim_output_dir = sim_output_path
        self.time_step = time_step
        self.object_data = self.get_object_data()

    def get_object_data(self):  # pylint: disable=unused-argument
        """Get aerodynamics data over an object inside the WindTunnel.
        
        Current Support - OpenFOAM
        """

        reading_file = os.path.join(self.sim_output_dir, "foam.foam")

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
        
        Returns:
            A MeshData object that allow to manipulate the data over a mesh
            and to render it.
        """

        property_notation = _get_physical_property_notation(physical_property)
        physical_field = MeshData(self.object_data, property_notation)

        return physical_field

    def get_streamlines(self, physical_property: str = "pressure"):
        """Get streamlines over the object in the WindTunnel."""

        streamlines_path = os.path.join(self.sim_output_dir, "postProcessing",
                                        "sets", "streamLines",
                                        str(self.time_step))

        property_notation = _get_physical_property_notation(physical_property)
        streamlines_file = "track0" + property_notation + ".vtk"

        streamlines_mesh = pv.read(
            os.path.join(streamlines_path, streamlines_file))

        return streamlines_mesh

    def get_flow_plane(self, physical_property: str = "pressure"):
        """Get flow properties in a plane of the domain in WindTunnel."""

        cutting_plane_path = os.path.join(self.sim_output_dir, "postProcessing",
                                          "cuttingPlane", str(self.time_step))

        # Obtain notation for the physical property for the simulator.
        property_notation = _get_physical_property_notation(physical_property)
        cutting_plane_file = property_notation + "_yNormal.vtk"

        cutting_plane_mesh = pv.read(
            os.path.join(cutting_plane_path, cutting_plane_file))

        return cutting_plane_mesh

    def render_flow(self,
                    flow_property_mesh,
                    physical_property: str = "pressure"):
        """Render flow property over the object in the WindTunnel."""

        # Obtain notation for the physical property for the simulator.
        property_notation = _get_physical_property_notation(physical_property)

        plotter = pv.Plotter()
        plotter.add_mesh(self.object_data)
        plotter.add_mesh(flow_property_mesh, scalars=property_notation)
        plotter.show()
        plotter.close()


def _get_physical_property_notation(physical_property: str):
    "Get OpenFOAM names for physical properties."

    if physical_property.lower() == "pressure":
        openfoam_property = "p"
    elif physical_property.lower() == "velocity":
        openfoam_property = "U"
    else:
        raise ValueError(
            f"Physical property `{physical_property}` not available "
            "or supported")

    return openfoam_property
