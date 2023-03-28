"""Utilities for creating mesh files for the fluid tank scenario."""

from . import gmsh_utils

from inductiva.fluids.shapes import Cube
from inductiva.fluids.shapes import Cylinder


def create_tank_mesh_file(shape, outlet, path: str):
    """Creates a mesh file for the tank.

    The tank is composed of two blocks:
    - a main (cylindrical/cubic) block representing the tank itself;
    - an optional smaller (cylindrical/cubic) block representing a fluid outlet.
      When present, the top base of this block connects with the bottom base of
      the tank, such that fluid flows freely from the tank to the outlet. The
      bottom base of the outlet is also open, such that flow exits the outlet.

    Both blocks are assumed to have their main axes aligned with the z
    axis.

    Args:
        shape: Shape of the tank.
        outlet: Shape of the outlet. If `None`, no outlet is present.
        path: Path of the file to be created.
    """

    with gmsh_utils.gmshAPIWrapper():
        tank_base_hole_loops = []

        if outlet is not None:

            # Add a circle arc/rectangle loop representing the top base of the
            # outlet. An arc/loop is used instead of a circle/rectangle because
            # this face is not filled, i.e. it is not a surface.
            if isinstance(outlet.shape, Cylinder):
                p_top_outlet, c_top_outlet, l_top_outlet = \
                    gmsh_utils.add_circle_arc(
                        x=outlet.shape.position[0],
                        y=outlet.shape.position[1],
                        z=outlet.shape.position[2] + outlet.shape.height,
                        r=outlet.shape.radius,
                    )
            elif isinstance(outlet.shape, Cube):
                p_top_outlet, c_top_outlet, l_top_outlet = \
                    gmsh_utils.add_z_rectangle_loop(
                        x=outlet.shape.position[0],
                        y=outlet.shape.position[1],
                        z=outlet.shape.position[2] + outlet.shape.dimensions[2],
                        lx=outlet.shape.dimensions[0],
                        ly=outlet.shape.dimensions[1],
                    )

            # Add a circle arc/rectangle loop representing the bottom base of
            # the outlet.
            if isinstance(outlet.shape, Cylinder):
                p_bottom_outlet, c_bottom_outlet, _ = gmsh_utils.add_circle_arc(
                    x=outlet.shape.position[0],
                    y=outlet.shape.position[1],
                    z=outlet.shape.position[2],
                    r=outlet.shape.radius,
                )
            elif isinstance(outlet.shape, Cube):
                p_bottom_outlet, c_bottom_outlet, _ = \
                    gmsh_utils.add_z_rectangle_loop(
                        x=outlet.shape.position[0],
                        y=outlet.shape.position[1],
                        z=outlet.shape.position[2],
                        lx=outlet.shape.dimensions[0],
                        ly=outlet.shape.dimensions[1],
                    )

            # Add the walls of the outlet (cylindrical/cubic) block.
            gmsh_utils.add_cylinder_walls(p_bottom_outlet, c_bottom_outlet,
                                          p_top_outlet, c_top_outlet)

            # Add the loop representing the top base of the outlet to the list
            # of loops representing holes in the bottom base of the tank
            # cylinder.
            tank_base_hole_loops.append(l_top_outlet)

        # Add the top and bottom bases of the tank block, setting the loop
        # representing the top base of the outlet as a hole.
        if isinstance(shape, Cylinder):
            p_top, c_top, _, _ = gmsh_utils.add_circle(
                x=shape.position[0],
                y=shape.position[1],
                z=shape.position[2] + shape.height,
                r=shape.radius,
                hole_loops=[],
            )
            p_bottom, c_bottom, _, _ = gmsh_utils.add_circle(
                x=shape.position[0],
                y=shape.position[1],
                z=shape.position[2],
                r=shape.radius,
                hole_loops=tank_base_hole_loops,
            )

        elif isinstance(shape, Cube):
            p_top, c_top, _, _ = gmsh_utils.add_z_rectangle(
                x=shape.position[0],
                y=shape.position[1],
                z=shape.position[2] + shape.dimensions[2],
                lx=shape.dimensions[0],
                ly=shape.dimensions[1],
                hole_loops=[],
            )
            p_bottom, c_bottom, _, _ = gmsh_utils.add_z_rectangle(
                x=shape.position[0],
                y=shape.position[1],
                z=shape.position[2],
                lx=shape.dimensions[0],
                ly=shape.dimensions[1],
                hole_loops=tank_base_hole_loops,
            )

        # Add the walls of the tank (cylindrical/cubic) block.
        gmsh_utils.add_cylinder_walls(p_bottom, c_bottom, p_top, c_top)

    # Convert the msh file generated by gmsh to obj format.
    gmsh_utils.convert_msh_to_obj_file(path)


def create_tank_fluid_mesh_file(shape, fluid_level, margin, path: str):
    """Creates a mesh file for the fluid.

    The fluid is represented by a block with the same shape as the tank, but
    with a smaller height.

    Args:
        shape: Shape of the tank.
        fluid_level: Height of the fluid.
        margin: Margin to be added to the fluid block.
        path: Path of the file to be created.
    """

    if isinstance(shape, Cube):
        with gmsh_utils.gmshAPIWrapper():
            gmsh_utils.add_box(
                shape.position[0] + margin,
                shape.position[1] + margin,
                shape.position[2] + margin,
                shape.dimensions[0] - 2 * margin,
                shape.dimensions[1] - 2 * margin,
                fluid_level - 2 * margin,
            )

    elif isinstance(shape, Cylinder):
        with gmsh_utils.gmshAPIWrapper():
            gmsh_utils.add_cylinder(
                shape.position[0],
                shape.position[1],
                shape.position[2] + margin,
                shape.radius - margin,
                fluid_level - margin,
            )
    else:
        raise ValueError(f"Invalid fluid shape `{shape}`.")

    # Convert the msh file generated by gmsh to obj format.
    gmsh_utils.convert_msh_to_obj_file(path)
