# Visualizing DualSPHysics Simulations with Blender

Blender is a powerful and versatile open-source 3D creation tool used by
professionals and hobbyists alike to craft stunning visualizations, animations,
and models. With its advanced rendering engine, support for complex shaders, and
a large ecosystem of plugins and tools, Blender is well-suited for creating
high-quality graphics and animations for everything from films and games to
scientific visualizations.

However, using Blender to visualize simulation results from DualSPHysics is not
straightforward. The simulation output, which typically consists of particle
data stored in vtk files, is not immediately compatible with Blender. To create
visualizations, users must perform post-processing steps—such as converting the
data into mesh representations or volumetric formats that Blender can understand.
This often involves using scripts or external tools to bridge the gap between
the raw simulation output and Blender's rendering system.


## Converting VTK Files to Blender-Compatible Formats

To visualize DualSPHysics simulations in Blender, the raw particle data
(typically stored in `.vtk` files) must first be converted into a mesh format
that Blender understands, such as `.obj`. This process involves several steps:
- Reading the particle data
- Constructing a volumetric density field
- Extracting a surface mesh from that field
- Exporting the mesh to an OBJ file.

Below is a Python script that automates this process. It reads `.vtk` files
containing particle positions, creates a 3D density field using a Gaussian kernel,
extracts an isosurface using the marching cubes algorithm, and saves the result
as a `.obj` file that can be imported directly into Blender.

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

Here's a clearer, more polished version of your tutorial. I've improved the flow, readability, and tone while keeping the structure and technical details intact:


## Rendering the Mesh in Blender

After generating your `.obj` files, you can import them into Blender by simply
dragging them into the scene. This will create one object for each frame of your
animation, resulting in many `.obj` files in your scene. While this isn't the
standard Blender workflow for animations, it's workable — we just need to
prepare the scene properly.

> **Note**: When you open Blender, feel free to delete the default cube from the
Scene Collection.

### 1. Set the Render Engine and Material

Since we'll be using ray tracing, first switch the render engine to **Cycles**:

* Go to the **Render** tab in the Properties panel.
* Set the **Render Engine** to `Cycles`.
* Set the **Device** to `GPU Compute` (if available).
<div align="center">
![Select Cycles for renderer](../../_static/render_engine.png)
<p align="center"><em>Figure 1: Select Cycles for renderer</em></p>
</div>
Next, assign a material to your objects:

1. Select the first object in the scene — for example, `PartFluid_0000`.
2. Go to the **Material** tab in the Properties panel.
3. Click **New** to create a new material. Rename it to something like `Water`.
4. Set the material properties:

   * **Base Color**: Light blue
   * **Roughness**: `0.1`
   * **IOR** (Index of Refraction): `1.33`
   * **Transmission**: `1.0` (for full transparency)
<div align="center">
![Creating the material](../../_static/material.png)
<p align="center"><em>Figure 2: Creating the material</em></p>
</div>
To apply this material to all the `.obj` objects:

1. Select `PartFluid_0000`, then Shift-click to select the last object (e.g., `PartFluid_0200`) to select all.
2. Press `Ctrl + L` and choose **Link Materials**.
<div align="center">
![Link the material to all objects](../../_static/material_all.png)
<p align="center"><em>Figure 3: Link the material to all objects</em></p>
</div>
Now all your fluid objects share the same water material. To preview the look,
change the viewport shading (top-right corner) to **Material Preview**.
<div align="center">
![Change the viewport](../../_static/viewport.png)
<p align="center"><em>Figure 4: Change the viewport</em></p>
</div>
### 2. Position the Camera

To set up the camera:

1. Select the camera in the **Scene Collection**.
2. In the **Object** tab of the Properties panel, set:

   * **Location**: `X: 200`, `Y: -120`, `Z: 650`
   * **Rotation**: `X: 0.0`, `Y: 0.0`, `Z: 180.0`
<div align="center">
![Position the camera](../../_static/blender_camera.png)
<p align="center"><em>Figure 5: Position the camera</em></p>
</div>
Additionally, in the **Camera Data** tab, set the **End** value of the lens to `1000`.
<div align="center">
![Set camera End](../../_static/camera_end.png)
<p align="center"><em>Figure 6: Set camera End</em></p>
</div>
### 3. Set Up the Animation

You might notice something odd: all frames are visible at once. That's because
we've imported each animation frame as a separate object. What we want is to
show only the relevant object for each frame.

To fix this, we'll write a small Blender script:

1. Go to the **Scripting** workspace (top bar).
2. Create a new script and paste the following:
<div align="center">
![Create a blender script](../../_static/scripting.png)
<p align="center"><em>Figure 7: Create a blender script</em></p>
</div>
```python
import bpy

# Select all fluid mesh objects
imported_objs = [obj for obj in bpy.data.objects if obj.type == 'MESH' and obj.name.startswith("PartFluid")]
imported_objs.sort(key=lambda o: o.name)  # Sort by name to ensure frame order

# Hide all objects initially
for obj in imported_objs:
    obj.hide_viewport = True
    obj.hide_render = True
    obj.keyframe_insert(data_path="hide_viewport", frame=1)
    obj.keyframe_insert(data_path="hide_render", frame=1)

# Animate: show one object per frame
for i, obj in enumerate(imported_objs):
    frame = i + 1

    obj.hide_viewport = False
    obj.hide_render = False
    obj.keyframe_insert(data_path="hide_viewport", frame=frame)
    obj.keyframe_insert(data_path="hide_render", frame=frame)

    obj.hide_viewport = True
    obj.hide_render = True
    obj.keyframe_insert(data_path="hide_viewport", frame=frame + 1)
    obj.keyframe_insert(data_path="hide_render", frame=frame + 1)
```

This script will:

* Hide all objects by default.
* Animate them so that each one is visible only on its corresponding frame.

### 4. Render the Animation

To export the animation:

1. Go to the **Output** tab in the Properties panel.
2. Set:

   * **Output directory**
   * **File format**
   * **End Frame** to `200`
<div align="center">
![Set the blender output](../../_static/output.png)
<p align="center"><em>Figure 8: Set the blender output</em></p>
</div>
Then, in the top menu, go to **Render > Render Animation**. Blender will render
each frame using ray tracing — this might take some time depending on your hardware.

### Wrapping Up

There's more you can explore, like adding an HDRI environment to enhance
lighting and realism, but we'll save that for another tutorial. For now, feel
free to tweak the material and camera settings to achieve your desired visual
effect.

**Happy rendering!**
