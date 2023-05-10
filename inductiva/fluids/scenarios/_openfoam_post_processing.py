"""Post process OpenFOAM simulation outputs.

Visualization Features:
- Scalar fields over object.
- Streamlines with scalar field coloring
- Cutting plane with scalar field coloring.
"""
import os
import pyvista as pv

from inductiva.types import Path


class OpenFOAMSimulationOutput:
    """Post process WindTunnel simulation outputs."""

    def __init__(self, sim_output_path: Path):
        """Initializes a `OpenFOAMSimulationOutput` object.

        Args:
            sim_output_path: Path to simulation output files.
        
        Attributes:
            sim_output_dir: path to the simulation directory
            object_path: path to the object insert in the WindTunnel
            post_processing_path: path to the post-processing data.
        """

        self.sim_output_dir = sim_output_path
        self.object_path = os.path.join(
            sim_output_path, "constant",
            "triSurface", "object.obj")
        self.post_processing_path = os.path.join(
            sim_output_path, "postProcessing")

    def render_scalar_field(self,
                            time_step: int,
                            scalars: str ="p"):
        """Render a scalar field over an object from OpenFOAM data.
        
        Args:
            time_step: The time frame for which we want to plot
                the scalar field over the object.
            scalars: The scalar field to plot. Meaningfull ones are
                the pressure "p" and velocity magnitude "U".
        """

        reading_file = os.path.join(self.sim_output_dir,
                                    "foam.foam")

        # Create the reading file for the WindTunnel scenarios
        open(reading_file, "w", encoding="utf-8").close()

        # Initialize reader and define reading time-step
        reader = pv.POpenFOAMReader(reading_file)
        reader.set_active_time_value(time_step)

        # Read mesh and mesh properties
        mesh = reader.read()
        object_mesh = mesh[1][5:]

        object_mesh.plot(scalars=scalars)

    def render_streamlines(self,
                           time_step: int,
                           scalars: str = "p"):
        """Render streamlines flowing by object with scalar color map.

        The streamlines are immediately obtained as outputs of 
        
        Args:
            time_step: The time frame for which we want to plot
                the scalar field over the object.
            scalars: The scalar field used to color the streamlines.
                Meaningfull ones are the pressure "p"
                and velocity magnitude "U".
        """

        object_mesh = pv.read(self.object_path)

        streamlines_path = os.path.join(
            self.post_processing_path, "sets", "streamLines",
            str(time_step))

        streamlines_file = "track0_" + scalars + ".vtk"
        streamlines_mesh = pv.read(os.path.join(streamlines_path,
                                           streamlines_file))

        plot = pv.Plotter()
        plot.add_mesh(streamlines_mesh, scalars=scalars)
        plot.add_mesh(object_mesh)

        plot.show()

    def render_cutting_plane(self,
                             time_step: int,
                             scalars: str = "p"):
        """Render the cutting plane obtain from the simulation.
        
        The cutting plane contains the scalar field properties over
        a plane over the domain. It is easier to see some of the
        properties of the flow.

        For now the only cutting plane being computed is orthogonal to
        the y-axis.

        Args:
            time_step: The time frame for which we want to plot
                the scalar field over the object.
            scalars: The scalar field used to color the streamlines.
                Meaningfull ones are the pressure "p"
                and velocity magnitude "U".
        """
        object_mesh = pv.read(self.object_path)

        cutting_plane_path = os.path.join(
            self.post_processing_path, "cuttingPlane",
            str(time_step))

        cutting_plane_file = scalars + "_yNormal.vtk"
        cutting_plane_mesh = pv.read(os.path.join(
            cutting_plane_path, cutting_plane_file))

        plot = pv.Plotter()
        plot.add_mesh(cutting_plane_mesh, scalars=scalars)
        plot.add_mesh(object_mesh)

        plot.show()
