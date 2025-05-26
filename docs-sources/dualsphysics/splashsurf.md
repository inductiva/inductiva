# Splashsurf

(Splashsurf)[https://github.com/InteractiveComputerGraphics/splashsurf] is an
open-source Rust-based tool designed for reconstructing surface
meshes from Smoothed Particle Hydrodynamics (SPH) simulation data.

We have integrated **splashsurf** as a post-processing tool that enables users
to convert VTK files containing SPH particle data into high-quality surface
meshes. These meshes will be saved as `.obj` files, making it easy for users to
import and render their simulation results in visualization tools such as Blender.
This integration streamlines the workflow from simulation to visualization,
providing a powerful bridge between numerical data and rendering.

## How to use splashsurf

We aim to make the process of using splashsurf as straightforward as possible.
This means that you only need a few extra arguments when running your simulation
to enable the splashsurf post-processing step. Here's how you can do it:

```diff
task = dualsphysics.run( \
    input_dir=input_dir,
    shell_script="xCaseDambreak3D_FSI_linux64_GPU.sh",
+    vtk_to_obj=True,
+    vtk_dir="CaseDambreak3D_FSI_out/particles/",
+    out_dir="CaseDambreak3D_FSI_out/particles/",
+    prefix="PartFluid_",
+    particle_radius=0.002,
+    smoothing_length=2,
+    cube_size=1,
+    surface_threshold=0.6,
    on=cloud_machine)
```


### Splashsurf Post-Processing Options

* **`vtk_to_obj`** *(bool)*:
  Enables the conversion of VTK particle data files into surface meshes using the marching cubes algorithm.

* **`vtk_dir`** *(str, mandatory if `vtk_to_obj=True`)*:
  Path to the directory containing the input VTK files.

* **`out_dir`** *(str, default: same as `vtk_dir`)*:
  Directory where the generated `.obj` mesh files will be saved. If not provided,
  saves the `.obj` files in the same folder as the `.vtk`files

* **`prefix`** *(str, default: `"PartFluid_"`)*:
  Prefix of the VTK files to be converted. Only files starting with this prefix
  will be processed. This will iterate over `PartFluid_1.vtk`, `PartFluid_2.vtk`
  etc.

* **`particle_radius`** *(float, mandatory if `vtk_to_obj=True`)*:
  Radius of each SPH particle in the input data.

* **`smoothing_length`** *(float, default: `2.0`)*:
  Smoothing length used in the SPH kernel, defined in multiples of the particle radius. The kernel’s compact support radius is twice this value.

* **`cube_size`** *(float, default: `1.0`)*:
  Edge length of the marching cubes grid cells, also in multiples of the particle radius. This defines the resolution of the reconstructed surface.

* **`surface_threshold`** *(float, default: `0.6`)*:
  Iso-surface threshold for the fluid density field, expressed in multiples of the SPH rest density. Determines the surface level used for mesh extraction.


Happy simulations!
