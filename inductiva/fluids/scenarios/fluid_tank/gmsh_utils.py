"""Utilities to parametrically generate meshes with gmsh.

TODO: Move this elsewhere in this repo to be reused by other high-level
scenarios? Or, further in the future, to API endpoints?
"""

import contextlib
import os

import gmsh
import meshio

TMP_MSH_FILE_PATH = "gmsh_generated.msh"
MESH_SIZE = 0.05
GMSH_LOG_LEVEL = 2


class gmshAPIWrapper(contextlib.AbstractContextManager):  # pylint: disable=invalid-name
    """Wrapper to interact with the gmsh API.

    Applies common operations with the gmsh API before and after use
    case-specific primitives.
    """

    def __init__(self,
                 mesh_size: float = MESH_SIZE,
                 msh_file_path: str = TMP_MSH_FILE_PATH,
                 log_level: int = GMSH_LOG_LEVEL):

        self.mesh_size = mesh_size
        self.msh_file_path = msh_file_path
        self.log_level = log_level

    def __enter__(self):
        """Common operations applied before use case-specific primitives."""

        # Initialize the gmsh API.
        gmsh.initialize()

        # Set the log level.
        gmsh.option.setNumber("General.Verbosity", self.log_level)

        # Set the mesh size.
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.mesh_size)

    def __exit__(self, *_):
        """Common operations applied after use case-specific primitives."""

        # Synchronize the representations (or kernels) available in gmsh:
        # - the built-in representation, used in a bottom-up manner (defining
        #   points first, then curves, surfaces and volumes).
        # - the OpenCascade representation, used as the built-in representation
        #   or in a top-down constructive solid geometry fashion (defining
        #   solids on which boolean operations are performed).
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        # Generate the mesh elements up to dimension 2, i.e. points, curves and
        # surfaces.
        gmsh.model.mesh.generate(2)

        # Write mesh to a file.
        gmsh.write(self.msh_file_path)

        # Close the gmsh API.
        gmsh.finalize()


def add_x_rectangle_loop(x, y, z, ly, lz):
    p1 = gmsh.model.geo.addPoint(x, y, z, 0)
    p2 = gmsh.model.geo.addPoint(x, y + ly, z, 0)
    p3 = gmsh.model.geo.addPoint(x, y + ly, z + lz, 0)
    p4 = gmsh.model.geo.addPoint(x, y, z + lz, 0)

    c1 = gmsh.model.geo.addLine(p1, p2)
    c2 = gmsh.model.geo.addLine(p2, p3)
    c3 = gmsh.model.geo.addLine(p3, p4)
    c4 = gmsh.model.geo.addLine(p4, p1)

    l = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])

    return [p1, p2, p3, p4], [c1, c2, c3, c4], l


def add_y_rectangle_loop(x, y, z, lx, lz):
    p1 = gmsh.model.geo.addPoint(x, y, z, 0)
    p2 = gmsh.model.geo.addPoint(x + lx, y, z, 0)
    p3 = gmsh.model.geo.addPoint(x + lx, y, z + lz, 0)
    p4 = gmsh.model.geo.addPoint(x, y, z + lz, 0)

    c1 = gmsh.model.geo.addLine(p1, p2)
    c2 = gmsh.model.geo.addLine(p2, p3)
    c3 = gmsh.model.geo.addLine(p3, p4)
    c4 = gmsh.model.geo.addLine(p4, p1)

    l = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])

    return [p1, p2, p3, p4], [c1, c2, c3, c4], l


def add_z_rectangle_loop(x, y, z, lx, ly):
    p1 = gmsh.model.geo.addPoint(x, y, z, 0)
    p2 = gmsh.model.geo.addPoint(x + lx, y, z, 0)
    p3 = gmsh.model.geo.addPoint(x + lx, y + ly, z, 0)
    p4 = gmsh.model.geo.addPoint(x, y + ly, z, 0)

    c1 = gmsh.model.geo.addLine(p1, p2)
    c2 = gmsh.model.geo.addLine(p2, p3)
    c3 = gmsh.model.geo.addLine(p3, p4)
    c4 = gmsh.model.geo.addLine(p4, p1)

    l = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])

    return [p1, p2, p3, p4], [c1, c2, c3, c4], l


def add_x_rectangle(x, y, z, ly, lz, hole_loops):

    # Add the corresponding circle arc.
    p, c, l = add_x_rectangle_loop(x, y, z, ly, lz)

    # Add hole loops to the list of loops defining the surface.
    ls = [l] + hole_loops

    # Add the surface.
    s = gmsh.model.geo.addPlaneSurface(ls)

    return p, c, l, s


def add_y_rectangle(x, y, z, lx, lz, hole_loops):

    # Add the corresponding circle arc.
    p, c, l = add_y_rectangle_loop(x, y, z, lx, lz)

    # Add hole loops to the list of loops defining the surface.
    ls = [l] + hole_loops

    # Add the surface.
    s = gmsh.model.geo.addPlaneSurface(ls)

    return p, c, l, s


def add_z_rectangle(x, y, z, lx, ly, hole_loops):

    # Add the corresponding circle arc.
    p, c, l = add_z_rectangle_loop(x, y, z, lx, ly)

    # Add hole loops to the list of loops defining the surface.
    ls = [l] + hole_loops

    # Add the surface.
    s = gmsh.model.geo.addPlaneSurface(ls)

    return p, c, l, s


def add_circle_arc(x, y, z, r):
    """Adds a circle arc to an open gmsh model.

    The circle arc is centered at [`x`, `y`, `z`], has radius `r` and is assumed
    to have a normal parallel to the z axis.
    """

    # First define the points required to draw the circle arc, i.e. the center
    # point (p1) and two pairs of diametrically opposed points (p2-p5).
    p1 = gmsh.model.geo.addPoint(x, y, z, 0)
    p2 = gmsh.model.geo.addPoint(x + r, y, z, 0)
    p3 = gmsh.model.geo.addPoint(x, y + r, z, 0)
    p4 = gmsh.model.geo.addPoint(x - r, y, z, 0)
    p5 = gmsh.model.geo.addPoint(x, y - r, z, 0)

    # Connect the multiple points into curves, e.g. c1 connects p2 to p3 with a
    # circle arc centered in p1.
    c1 = gmsh.model.geo.addCircleArc(p2, p1, p3)
    c2 = gmsh.model.geo.addCircleArc(p3, p1, p4)
    c3 = gmsh.model.geo.addCircleArc(p4, p1, p5)
    c4 = gmsh.model.geo.addCircleArc(p5, p1, p2)

    # Merge the multiple curves into a loop.
    l = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])

    return [p2, p3, p4, p5, p1], [c1, c2, c3, c4], l


def add_circle(x, y, z, r, hole_loops):
    """Adds a circle to an open gmsh model.

    The circle is centered at [`x`, `y`, `z`], has radius `r` and is assumed to
    have a normal parallel to the z axis.

    A circle differs from a circle arc in that it has its interior filled with a
    (plane) surface. This surface can have holes represented by the loops in
    `hole_loops`.
    """

    # Add the corresponding circle arc.
    p, c, l = add_circle_arc(x, y, z, r)

    # Add hole loops to the list of loops defining the surface.
    ls = [l] + hole_loops

    # Add the surface.
    s = gmsh.model.geo.addPlaneSurface(ls)

    return p, c, l, s


def add_cylinder_walls(p_base, c_base, p_top, c_top):
    """Adds cylinder walls (i.e. the curved surface) to an open gmsh model.

    The cylinder walls are built by merging the points and curves defining the
    bottom and top bases of the cylinder.

    More specifically, they are built in four steps, each focused on a curve
    defining a quarter of the base circle arcs.

    In each step, a line is created to connect a point in the bottom base with
    the corresponding point in the top base. Then, the sub-arc of the top base
    starting at that point is connected with the line. Another line is then
    added to connect the end point of the sub-arc in the top base to the
    corresponding point in the bottom base. Finally, the missing sub-arc in the
    bottom base is added (reversed relative to its definition, to ensure that
    the final set of curves defines a closed loop).
    """

    num_curves = len(c_base)

    # Loop through the steps, which number matches the number of curves in the
    # bases of the cylinder.
    for idx in range(num_curves):

        # Define the start and end indices.
        idx_start = idx
        idx_end = (idx + 1) % num_curves

        # Get the list of curves defining a quarter of the cylinder curved
        # surface:
        # 1. a line connecting the bottom and top bases;
        # 2. a sub-arc in the top base;
        # 3. a line connecting the top and bottom bases;
        # 4. a sub-arc in the bottom base (reversed relative to its definition,
        #    to ensure that the final set of curves defines a closed loop)-
        curves = [
            gmsh.model.geo.addLine(p_base[idx_start], p_top[idx_start]),
            c_top[idx],
            gmsh.model.geo.addLine(p_top[idx_end], p_base[idx_end]),
            -c_base[idx],
        ]

        loop = gmsh.model.geo.addCurveLoop(curves)
        gmsh.model.geo.addSurfaceFilling([loop])


def add_box(x, y, z, lx, ly, lz):
    gmsh.model.occ.addBox(
        x,
        y,
        z,
        lx,
        ly,
        lz,
    )


def add_cylinder(x, y, z, radius, height):
    gmsh.model.occ.addCylinder(
        x,
        y,
        z,
        0,
        0,
        height,
        radius,
    )


def convert_msh_to_obj_file(obj_file_path: str,
                            msh_file_path: str = TMP_MSH_FILE_PATH,
                            remove_msh_file: bool = True) -> None:
    """Converts msh to obj files using meshio."""

    # Read mesh using meshio.
    mesh = meshio.read(msh_file_path, file_format="gmsh")

    # Filter cell blocks to keep only those representing surfaces.
    cells = [block for block in mesh.cells if block.type == "triangle"]

    # Write mesh to obj file.
    meshio.write_points_cells(obj_file_path, mesh.points, cells)

    if remove_msh_file:
        os.remove(msh_file_path)
