# Converting VTK Files to Blender-Compatible Formats

To visualize DualSPHysics simulations in Blender, the raw particle data
(typically stored in `.vtk` files) must first be converted into a mesh format
that Blender understands, such as `.obj`. This process can be easily done with
Inductiva.

Let's use the [3D dam break](../../multiple_gpus.md) tutorial as an example. In
this simulation we used 

```python
import vtk
import numpy as np
from numba import njit, prange
from vtk.util.numpy_support import vtk_to_numpy
from scipy.ndimage import gaussian_filter
from skimage import measure
import os

def vtk_to_points(vtk_file):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtk_file)
    reader.Update()
    data = reader.GetOutput()
    points = vtk_to_numpy(data.GetPoints().GetData())
    return points

@njit(parallel=True)
def create_density_field(points, grid_size=500, radius=0.005):
    field = np.zeros((grid_size, grid_size, grid_size), dtype=np.float32)
    spacing = 1.0 / grid_size
    r2 = radius * radius
    sigma2 = (radius / 2) ** 2
    n_points = points.shape[0]

    for p_idx in prange(n_points):
        px, py, pz = points[p_idx]
        i, j, k = int(px / spacing), int(py / spacing), int(pz / spacing)
        for dx in range(-2, 3):
            for dy in range(-2, 3):
                for dz in range(-2, 3):
                    ni, nj, nk = i + dx, j + dy, k + dz
                    if 0 <= ni < grid_size and 0 <= nj < grid_size and 0 <= nk < grid_size:
                        gx = ni * spacing
                        gy = nj * spacing
                        gz = nk * spacing
                        dist2 = (gx - px)**2 + (gy - py)**2 + (gz - pz)**2
                        if dist2 < r2:
                            field[ni, nj, nk] += np.exp(-dist2 / (2 * sigma2))
    return field

def mesh_from_field(field, iso=0.5):
    verts, faces, _, _ = measure.marching_cubes(field, level=iso)
    return verts, faces

def save_obj(filename, verts, faces):
    with open(filename, 'w') as f:
        for v in verts:
            f.write(f"v {v[0]} {v[1]} {v[2]}\n")
        for face in faces:
            f.write(f"f {int(face[0]+1)} {int(face[1]+1)} {int(face[2]+1)}\n")

def convert_vtk_dir_to_meshes(vtk_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for file in sorted(os.listdir(vtk_dir)):
        print(f"Checking file {file}...")
        if file.endswith('.vtk') and file.startswith("PartFluid"):
            print("It's vtk and PartFluid ...")
            path = os.path.join(vtk_dir, file)
            print("Starting vtk to points...")
            points = vtk_to_points(path)
            print("Starting create density field...")
            field = create_density_field(points)
            print("Starting mesh from field...")
            verts, faces = mesh_from_field(field)
            frame_id = os.path.splitext(file)[0]
            save_obj(os.path.join(out_dir, f"{frame_id}.obj"), verts, faces)
            print(f"Converted {file} → mesh")

# Example usage
convert_vtk_dir_to_meshes(
    "/Path/to/particles/folder",
    "/Path/to/store/the/results"
)
```

This script automates the conversion of `.vtk` files—typically named like
`PartFluid_000000.vtk`—into `.obj` mesh files, which can be easily imported and
rendered in Blender. Since Blender cannot directly visualize particle-based data,
this conversion step is essential for creating meaningful and visually appealing
representations of DualSPHysics simulations.

The core of the script involves converting particle positions into a volumetric
density field using a Gaussian kernel, then extracting a surface using the
marching cubes algorithm. The resulting mesh is saved as an `.obj` file, ready
for use in Blender.

Several parameters in the script allow you to control the quality, appearance,
and performance of the mesh generation process:

* **`grid_size`** - Sets the resolution of the 3D density grid.
  Higher values produce smoother and more detailed surfaces but increase memory usage and runtime. Lower values are faster but result in coarser meshes.

* **`radius`** - Defines how far each particle influences the density field.
  A smaller radius results in sharper, more detailed features, while a larger radius yields smoother, more continuous surfaces.

* **`iso`** - The threshold used by the marching cubes algorithm to define the surface.
  Higher values generate thinner, more compact meshes; lower values produce thicker and more voluminous ones.

By tweaking these parameters, you can tailor the output for different use
cases—ranging from quick previews to high-quality final renders.

> **Note:** Make sure to install the required dependencies before running the script (e.g., `pip install vtk scikit-image numba`). Also, if you have multiple sets of `.vtk` files from different simulation outputs, you may need to run the script separately for each group.

See how you can visualize the results in Blender in the next section [Rendering in Blender](render_in_blender.md).
