"""Class interface to manipulate CADObjects.

With the goal to manipulate the geometry and spatial properties
of the objects to be used inside our simulators, this class
provides a simple interface to read the objects from a native
CAD file.
"""
from typing import List, Literal, Union

import pyvista as pv


class CADObject:

    def __init__(self, filename):
        """"Initialize the CADObject class."""

        self.filename = filename
        self.data = pv.read(filename)

    def __repr__(self):
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
            inplace: If True, the object is translated in place.
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
            inplace: If True, the object is translated in place.
        """

        print("Translating object to the origin.")
        data = self.data.translate([-value for value in self.data.center])
        
        if inplace:
            self.data = data

        return self
        

    def to_ground(self, inplace: bool = False):
        """Translate the CADObject to the ground.

        Args:
            inplace: If True, the object is translated in place.
        """

        print("Translating object to the ground.")
        data = self.data.translate([0., 0., -self.data.center[2]])

        if inplace:
            self.data = data

        return self

    def rotate_around_vector(self,
                             vector: List[float],
                             center: List[float] = [0., 0., 0.],
                             angle: float = 0.,
                             inplace: bool = False):
        """Rotate the CADObject by an angle around a vector.
        
        Args:
            vector: 3D vector to rotate the object around.
            center: 3D center point of the rotation.
            angle: Angle in degrees to rotate the object.
            inplace: If True, the object rotation is saved in the data.
        """

        check_spatial_vector(vector, "Vector")
        check_spatial_vector(center, "Center")

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
            angle: Angle in degrees to rotate the object.
            inplace: If True, the object rotation is saved in the data.
        """

        axis_vector_dict = {
            "x": [1., 0., 0.],
            "y": [0., 1., 0.],
            "z": [0., 0., 1.]
        }

        if axis in ["x", "y", "z"]:
            vector = axis_vector_dict[axis]
        else:
            raise ValueError("Axis must be one of 'x', 'y', 'z'.")

        print(f"Rotating around the {axis}-axis by {angle} degrees.")
        data = self.data.rotate_vector(vector, angle, inplace=inplace)

        if inplace:
            self.data = data

        return self
    
    def scale(self, xyz: Union[float, List[float]], inplace: bool = False):
        """Scale the CADObject by a vector.
        
        Args:
            vector: float or 3D vector that define the scale factors,
                overall and along x, y, and z, respectively. If a float
                is passed the same uniform scale is applied to all axes.
            inplace: If True, the object is scaled in place.
        """

        if type(xyz) == float:
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

    def show(self,
             show_edges: bool = False,
             off_screen: bool = False,
             virtual_display: bool = False,
             object_color: str = "white",
             background_color: str = "grey",
             save_path: str = None
             ):
        """Show the CADObject data in a PyVista plot."""

        if virtual_display:
            off_screen = True
            pv.start_xvfb()

        plotter = pv.Plotter(off_screen=off_screen)
        plotter.background_color = background_color
        plotter.add_camera_orientation_widget()

        plotter.add_mesh(self.data, show_edges=show_edges, color=object_color)
        plotter.show(auto_close=True)

        if save_path is not None:
            plotter.screenshot(save_path, return_img=False)
            print(f"CADObject render was saved to {save_path}")

        plotter.close()


def check_spatial_vector(vector: List[float], vector_name: str):
    """Check if vector represents a 3D spatial vector."""

    if len(vector) != 3:
        raise ValueError(f"{vector_name} must be three dimensional.")

    for value in vector:
        if type(value) != float:
            raise ValueError(f"{vector_name} has an invalid value format. "
                             "Vector needs to contain "
                             "only float values.")
