"""Utils to create the geometric case."""

from typing import List, Optional

import gmsh
import json
import math


class RectangularPlate:
    """Rectangular plate.

    Attributes:
        plate_length: A float representing the plate length.
        plate_width: A float representing the plate width.
    """

    def __init__(self, plate_length: float, plate_width: float) -> None:
        """Initializes a RectangularPlate object."""
        self.plate_length = plate_length
        self.plate_width = plate_width


class CircularHole:
    """Circular hole."""

    def __init__(self, center_x: float, center_y: float, radius: float) -> None:
        """Initializes a CircularHole object.

        Attributes:
            hole_type: A string representing the hole type.
            center_x: A float representing the x coordinate of hole center.
            center_y: A float representing the y coordinate of hole center.
            radius: A float representing the hole radius.
        """
        self.hole_type = "circular"
        self.center_x = center_x
        self.center_y = center_y
        self.radius = radius


class RectangularHole:
    """Rectangular hole."""

    def __init__(self, center_x: float, center_y: float, size_x: float,
                 size_y: float, angle: float) -> None:
        """Initializes a RectangularHole object.

        Attributes:
            hole_type: A string representing the hole type.
            center_x: A float representing the x coordinate of hole center.
            center_y: A float representing the y coordinate of hole center.
            size_x: A float representing the hole size in x direction.
            size_y: A float representing the hole size in y direction.
            angle: A float representing the positive angle of rotation in
              degrees around the hole center.
        """
        self.hole_type = "rectangular"
        self.center_x = center_x
        self.center_y = center_y
        self.size_x = size_x
        self.size_y = size_y
        self.angle = angle


class EllipticalHole:
    """Elliptical hole."""

    def __init__(self, center_x: float, center_y: float, radius_x: float,
                 radius_y: float, angle: float) -> None:
        """Initializes a EllipticalHole object.

        Attributes:
            hole_type: A string representing the hole type.
            center_x: A float representing the x coordinate of hole center.
            center_y: A float representing the y coordinate of hole center.
            radius_x: A float representing the hole radius in x direction.
            radius_y: A float representing the hole radius in y direction.
            angle: A float representing the positive angle of rotation in
              degrees around the hole center.
        """
        self.hole_type = "elliptical"
        self.center_x = center_x
        self.center_y = center_y
        self.radius_x = radius_x
        self.radius_y = radius_y
        self.angle = angle


class GeometricCase:
    """Geometric case.

    The geometric case is characterized by a plate and a set of holes.

    Attributes:
        plate_object: An instance of one of the RectangularPlate class.
        list_holes_object: A list of instances of holes classes:
          CircularHole, RectangularHole, or EllipticalHole.
    """

    def __init__(
        self,
        plate_object: RectangularPlate,
        list_holes_objects: Optional[List[CircularHole or RectangularHole or
                                          EllipticalHole]] = None
    ) -> None:
        """Initializes a GeometricCase object."""
        self.plate_object = plate_object
        self.list_holes_objects = list_holes_objects

    def write_to_json(self, json_path: str) -> None:
        """Writes the geometric case to JSON file.

        Args:
            json_path: A string representing the JSON file path.
        """

        # Plate dictionary
        plate_dict = {"plate": vars(self.plate_object)}

        # Holes dictionary
        if self.list_holes_objects is not None:

            all_hole_properties = []
            for hole_object in self.list_holes_objects:
                hole_dictionary = vars(hole_object)
                all_hole_properties.append(hole_dictionary)

            holes_dict = {"holes": all_hole_properties}

            # Merge dictionaries: geometric case dictionary
            geom_case_dictionary = {**plate_dict, **holes_dict}
        else:
            geom_case_dictionary = {**plate_dict}

        # Write JSON file
        with open(json_path, "w", encoding="utf-8") as write_file:
            json.dump(geom_case_dictionary, write_file, indent=4)

    def _plate_to_opencascade_cad(self) -> tuple:
        """Converts plate object to OpenCASCADE CAD representation.

        Raises:
            ValueError: If there is no plate object attributed to the 
              geometric case.

        Returns:
          tuple: A tuple containing the plate gmsh object and the IDs of the
            curves that constitute the plate.
        """
        if self.plate_object is not None:

            # Rectangular plate
            plate_gmsh = gmsh.model.occ.addRectangle(
                x=0,
                y=0,
                z=0,
                dx=self.plate_object.plate_width,
                dy=self.plate_object.plate_length)

            gmsh.model.occ.synchronize()

            plate_curves_id = []

            for curve_id in gmsh.model.getBoundary([(2, plate_gmsh)]):
                plate_curves_id.append([curve_id[1]][0])

            return plate_gmsh, plate_curves_id

        else:
            raise ValueError("There is no plate object attributed to the "
                             "geometric case. Please check.")

    def _get_plate_metrics_for_mesh_generation(self) -> Optional[tuple]:
        """Gets the plate metrics required for building the mesh in gmsh.

        Metrics:
          - mesh_offset (float): Represents an offset for all the curves of the
            plate, defining a region around the boundaries. Within this region, 
            we have the ability to control the mesh elements size.
            For rectangular plates, the offset will be equal to half of the
            minimum size of the plate.
          - predefined_element_size (float): Represents the predefined element
            size, defined as 1/4 of the perimeter.

        Raises:
            ValueError: If there is no plate object attributed to the 
              geometric case.

        Returns:
          tuple: A tuple containing the mesh offset and the predefined element
            mesh size for the plate.
      """
        if self.plate_object is not None:

            # Rectangular plate

            plate_mesh_offset = min(self.plate_object.plate_width,
                                    self.plate_object.plate_length) / 2
            plate_predefined_element_size = (
                self.plate_object.plate_width * 2 +
                self.plate_object.plate_length * 2) / 4

            return plate_mesh_offset, plate_predefined_element_size

        else:
            raise ValueError("There is no plate object attributed to the "
                             "geometric case. Please check.")

    def _holes_to_opencascade_cad(self) -> Optional[tuple]:
        """Converts holes objects to OpenCASCADE CAD representation.

        Raises:
            ValueError: If the hole type is different from the available
              options.

        Returns:
            tuple: A tuple containing the OpenCASCADE CAD representation of the
              holes and the IDs of the curves of the holes.
        """

        if self.list_holes_objects:
            holes_gmsh = []
            holes_curve_id = []

            for hole_object in self.list_holes_objects:

                # Circular hole
                if hole_object.hole_type == "circular":
                    hole_gmsh = gmsh.model.occ.addDisk(xc=hole_object.center_x,
                                                       yc=hole_object.center_y,
                                                       zc=0,
                                                       rx=hole_object.radius,
                                                       ry=hole_object.radius)

                # Rectangular hole
                elif hole_object.hole_type == "rectangular":
                    hole_gmsh = gmsh.model.occ.addRectangle(
                        x=hole_object.center_x - hole_object.size_x,
                        y=hole_object.center_y - hole_object.size_y,
                        z=0,
                        dx=hole_object.size_x * 2,
                        dy=hole_object.size_y * 2)
                    gmsh.model.occ.rotate(dimTags=[(2, hole_gmsh)],
                                          x=hole_object.center_x,
                                          y=hole_object.center_y,
                                          z=0,
                                          ax=0,
                                          ay=0,
                                          az=1,
                                          angle=math.radians(hole_object.angle))

                # Elliptical hole
                elif hole_object.hole_type == "elliptical":

                    if hole_object.radius_x > hole_object.radius_y:
                        hole_gmsh = gmsh.model.occ.addDisk(
                            xc=hole_object.center_x,
                            yc=hole_object.center_y,
                            zc=0,
                            rx=hole_object.radius_x,
                            ry=hole_object.radius_y)
                        gmsh.model.occ.rotate(dimTags=[(2, hole_gmsh)],
                                              x=hole_object.center_x,
                                              y=hole_object.center_y,
                                              z=0,
                                              ax=0,
                                              ay=0,
                                              az=1,
                                              angle=math.radians(
                                                  hole_object.angle))

                    else:
                        hole_gmsh = gmsh.model.occ.addDisk(
                            xc=hole_object.center_x,
                            yc=hole_object.center_y,
                            zc=0,
                            rx=hole_object.radius_y,
                            ry=hole_object.radius_x)
                        gmsh.model.occ.rotate(
                            dimTags=[(2, hole_gmsh)],
                            x=hole_object.center_x,
                            y=hole_object.center_y,
                            z=0,
                            ax=0,
                            ay=0,
                            az=1,
                            angle=math.radians(hole_object.angle + 90))

                else:
                    raise ValueError("The hole type is different from "
                                     "circular, rectangular, or elliptical. "
                                     "Please check.")

                gmsh.model.occ.synchronize()

                hole_curve_id = []
                for curve_id in gmsh.model.getBoundary([(2, hole_gmsh)]):
                    hole_curve_id.append([curve_id[1]][0])

                holes_gmsh.append(hole_gmsh)
                holes_curve_id.append(hole_curve_id)

            return holes_gmsh, holes_curve_id

    def _get_holes_metrics_for_mesh_generation(self) -> Optional[tuple]:
        """Gets the holes metrics required for building the mesh in gmsh. 

        Metrics:
          - mesh_offset (float): Represents an offset for the curves of the
            holes, defining a region around the boundaries. Within this region, 
            we have the ability to control the mesh elements size.
            For rectangular plates, the offset will be equal to half of the
            minimum size of the plate. For circular holes, it will be the
            radius. For elliptical holes, it will be the minimum radius.
          - predefined_element_size (float): Represents the predefined element
            size, defined as 1/4 of the perimeter.

        Raises:
            ValueError: If the hole type is different from the available
              options.

        Returns:
            tuple: A tuple containing the mesh offset and the predefined element
              mesh size for all the holes.
        """

        if self.list_holes_objects:

            holes_mesh_offset = []
            holes_predefined_element_size = []

            for hole_object in self.list_holes_objects:

                # Circular hole
                if hole_object.hole_type == "circular":

                    hole_mesh_offset = hole_object.radius
                    hole_predefined_element_size = (2 * math.pi *
                                                    hole_object.radius / 4)

                # Rectangular hole
                elif hole_object.hole_type == "rectangular":

                    hole_mesh_offset = min(hole_object.size_x,
                                           hole_object.size_y) / 2
                    hole_predefined_element_size = (
                        self.plate_object.plate_width * 2 +
                        self.plate_object.plate_length * 2) / 4

                # Elliptical hole
                elif hole_object.hole_type == "elliptical":

                    hole_mesh_offset = min(hole_object.radius_x,
                                           hole_object.radius_y)
                    hole_predefined_element_size = math.pi * (
                        3 * (hole_object.radius_x + hole_object.radius_y) -
                        math.sqrt(
                            (3 * hole_object.radius_x + hole_object.radius_y) *
                            (hole_object.radius_x + 3 * hole_object.radius_y))
                    ) / 4

                else:
                    raise ValueError("The hole type is different from "
                                     "circular, rectangular, or elliptical. "
                                     "Please check.")

                holes_mesh_offset.append(hole_mesh_offset)
                holes_predefined_element_size.append(
                    hole_predefined_element_size)

            return holes_mesh_offset, holes_predefined_element_size

    def to_gmsh(self) -> tuple:
        """Generates the plate with holes in the OpenCASCADE CAD representation.

        The process of generating the plate with holes is divided into three
        steps:

        1. Converts plate object to OpenCASCADE CAD representation
        2. Converts holes objects to OpenCASCADE CAD representation
        3. Removes holes from the plate in the OpenCASCADE CAD representation
        4. Gets the holes and plate metrics required for building the mesh

        To remove the holes from the plate, the gmsh.model.occ.cut() function
        will be used.

        Returns:
          tuple: A tuple containing the IDs of the curves, the mesh offset and
            the predefined element mesh size for the plate and all the holes.
        """

        # 1. Converts the plate object to OpenCASCADE CAD representation
        plate_gmsh, plate_curves_id = self._plate_to_opencascade_cad()

        # 2. Converts the hole objects to OpenCASCADE CAD representation
        holes_gmsh, holes_curve_id = self._holes_to_opencascade_cad()

        # 3. Removes holes from the plate in the OpenCASCADE CAD representation
        for hole in holes_gmsh:
            gmsh.model.occ.cut([(2, plate_gmsh)], [(2, hole)])

        gmsh.model.occ.synchronize()

        # 4. Gets the holes and plate metrics required for building the mesh
        (plate_mesh_offset, plate_predefined_element_size
         ) = self._get_plate_metrics_for_mesh_generation()
        (holes_mesh_offset, holes_predefined_element_size
         ) = self._get_holes_metrics_for_mesh_generation()

        return (plate_curves_id, plate_mesh_offset,
                plate_predefined_element_size, holes_curve_id,
                holes_mesh_offset, holes_predefined_element_size)
