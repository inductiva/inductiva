"""Class interface to manipulate CADObjects.

With the goal to manipulate the geometry and spatial properties
of the objects to be used inside our simulators, this class
provides a simple interface to read the objects from a native
CAD file.
"""
from typing import List, Literal, Optional, Union

try:
    import pyvista as pv
except ImportError:
    pv = None

from inductiva.utils import optional_deps


class CADObject:
    """Class to manipulate CADObjects.
    
    Computer-aided design (CAD) objects are used to represent
    the geometry of the objects to be simulated. This class
    provides a simple interface to read the objects from a
    native CAD file and manipulate them in space.
    """

    @optional_deps.needs_common_extra_deps
    def __init__(self, filename):
        """"Initialize the CADObject class.
        
        Args:
            filename (str): The path to the CAD file
                to be loaded. 
                Supported formats: .stl, .ply, .vtk, .obj,
                .vtu, .ply, .msh and more. Check the link
                https://github.com/nschloe/meshio for more.
        """

        self.filename = filename
        self.data = pv.read(filename)

    def __repr__(self) -> str:
        """Print information about the CAD object.
        
        TODO: Add more information about the object.
        Existing regions, materials, etc...
        """

        print("\nCADObject details: ")
        print(self.data)
        print(f"\nMesh center: {self.data.center}\n")

        return f"CADObject filename: {self.filename}"

    def translate(self, vector: List[float], inplace: bool = False):
        """Translate the CADObject by a vector.

        Args:
            vector: The vector to translate the object.
            inplace: If True, the object is saved in the data.

        Returns:
            The CADObject instance after translating.
        """
        check_spatial_vector(vector, "Vector")

        print("Translate object by: ", vector)
        data = self.data.translate(vector, inplace=inplace)

        if inplace:
            self.data = data

        return self

    def to_origin(self, inplace: bool = False):
        """Translate the CADObject to the origin.

        Args:
            inplace: If True, the object is save in the data.

        Returns:
            The CADObject instance centered in the origin.
        """

        print("Translating object to the origin.")
        data = self.data.translate([-value for value in self.data.center])

        if inplace:
            self.data = data

        return self

    def rotate_around_vector(self,
                             vector: List[float],
                             center: Optional[List[float]] = None,
                             angle: float = 0.,
                             inplace: bool = False):
        """Rotate the CADObject by an angle around a vector.
        
        Args:
            vector: 3D vector to rotate the object around.
            center: 3D center point of the rotation.
                Default is None, which is converted to [0., 0., 0.]
            angle: Angle in degrees to rotate the object.
            inplace: If True, the object rotation is saved in the data.

        Returns:
            The CADObject instance after rotating around a vector.
        """

        check_spatial_vector(vector, "Vector")
        check_spatial_vector(center, "Center")

        if center is None:
            center = [0., 0., 0.]

        print(f"Rotating around {vector} at the center point {center}"
              f" by {angle} degrees")
        data = self.data.rotate_vector(vector, angle, center, inplace)

        if inplace:
            self.data = data

        return self

    def rotate_around_axis(self,
                           axis: Literal["x", "y", "z"],
                           angle: float = 0.,
                           inplace: bool = False):
        """Rotate the CADObject by an angle around an cartesian axis.
        
        Args:
            axis: Axis to rotate the object around. 
                Available options: "x", "y" and "z".
            angle: Angle in degrees to rotate the object.
            inplace: If True, the object rotation is saved in the data.

        Returns:
            The CADObject instance after rotating around an axis.
        """

        axis_vector_dict = {
            "x": [1., 0., 0.],
            "y": [0., 1., 0.],
            "z": [0., 0., 1.]
        }

        vector = axis_vector_dict[axis]
        print(f"Rotating around the {axis}-axis by {angle} degrees.")
        data = self.data.rotate_vector(vector, angle, inplace=inplace)

        if inplace:
            self.data = data

        return self

    def scale(self, xyz: Union[float, List[float]], inplace: bool = False):
        """Scale the CADObject by a vector.
        
        Args:
            xyz: float or 3D vector that define the scale factors,
                overall and along x, y, and z, respectively. If a float
                is passed the same uniform scale is applied to all axes.
            inplace: If True, the object is scaled in place.

        Returns:
            The CADObject instance after scaling.
        """

        if isinstance(xyz, float):
            vector = [xyz, xyz, xyz]
            print("Scaling object by a factor of ", xyz)
        else:
            check_spatial_vector(xyz, "Vector")
            print("Scaling object with scale factors ", xyz)
            vector = xyz

        data = self.data.scale(vector, inplace=inplace)

        if inplace:
            self.data = data

        return self

    def save(self, filename: str, binary: bool = False):
        """Save the current mesh to a file.
        
        Args:
            filename (str): The path to the file to be saved.
                Available formats: .stl, .ply. and .vtk.
            binary (bool): Whether to save the file in binary format.
                If False, the file is saved in ASCII format.
        """

        file_format = filename.split(".")[-1]

        if file_format not in ["stl", "ply", "vtk"]:
            raise ValueError("File format not supported. "
                             "Available formats: .stl, .ply, .vtk")

        self.data.save(filename, binary)

    def show(self,
             show_edges: bool = False,
             off_screen: bool = False,
             virtual_display: bool = False,
             object_color: str = "white",
             background_color: str = "grey",
             n_ticks: int = 2,
             save_path: str = None):
        """Show the CADObject data in a PyVista plot.
        
        Args:
        Args:
            show_edges (bool): Whether to show edges of the CAD object.
            off_screen (bool): Whether to render off-screen.
            virtual_display (bool): Whether to use a virtual display
                (requires xvfb).
            object_color (str): Color of the CAD object.
            background_color (str): Background color of the render.
            save_path (str, optional): If provided, the render will be
                saved to this path.
        """

        if virtual_display:
            off_screen = True
            pv.start_xvfb()

        plotter = pv.Plotter(off_screen=off_screen)
        plotter.background_color = background_color
        plotter.add_camera_orientation_widget()

        plotter.add_mesh(self.data, show_edges=show_edges, color=object_color)
        # Set axis on the back of the object and add a padding to the grid
        plotter.show_grid(location="outer",
                          padding=0.3,
                          n_xlabels=n_ticks,
                          n_ylabels=n_ticks,
                          n_zlabels=n_ticks,
                          font_size=9)

        plotter.show(auto_close=True)

        if save_path is not None:
            plotter.screenshot(save_path, return_img=False)
            print(f"CADObject render was saved to {save_path}")

        plotter.close()


def check_spatial_vector(vector: List[float], vector_name: str):
    """Check if vector represents a 3D spatial vector.
    
    Args:
        vector (List[float]): The list to be checked for 3D spatial 
            vector representation.
        vector_name (str): A descriptive name for the vector for error
            messages.

Raises:
        ValueError: If the input list does not have exactly three
            float values.
            If any element in the list is not a float."""

    if len(vector) != 3:
        raise ValueError(f"{vector_name} must be three dimensional.")

    for value in vector:
        if not isinstance(value, float):
            raise ValueError(f"{vector_name} has an invalid value format. "
                             "Vector needs to contain "
                             "only float values.")
