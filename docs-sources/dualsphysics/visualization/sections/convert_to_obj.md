# Converting VTK Files to Blender-Compatible Formats

To visualize DualSPHysics simulations in Blender, the raw particle data
(typically stored in `.vtk` files) must first be converted into a mesh format
that Blender understands, such as `.obj`. This process can be easily done with
Inductiva.

Let's use the [3D dam break](../../multiple_gpus.md) tutorial as an example. In
this simulation we run this simulation with the following parameters:

```python
# Run simulation
task = dualsphysics.run( \
    input_dir=input_dir,
    shell_script="xCaseDambreak3D_FSI_linux64_GPU.sh",
    vtk_to_obj=True,
    vtk_to_obj_vtk_dir="CaseDambreak3D_FSI_out/particles/",
    vtk_to_obj_vtk_prefix="PartFluid_",
    vtk_to_obj_particle_radius=0.002,
    vtk_to_obj_smoothing_length=2,
    vtk_to_obj_cube_size=1,
    on=cloud_machine)
```

This will use a tool called **splashsurf** to convert the VTK files into `.obj`.

## Splashsurf

(Splashsurf)[https://github.com/InteractiveComputerGraphics/splashsurf] is an
open-source Rust-based tool designed for reconstructing surface
meshes from Smoothed Particle Hydrodynamics (SPH) simulation data.

We have integrated **splashsurf** as a post-processing tool that enables users
to convert VTK files containing SPH particle data into high-quality surface
meshes. These meshes will be saved as `.obj` files, making it easy for users to
import and render their simulation results in visualization tools such as Blender.
This integration streamlines the workflow from simulation to visualization,
providing a powerful bridge between numerical data and rendering.

### Splashsurf Post-Processing Options

* **`vtk_to_obj`** *(bool)*:
  Enables the conversion of VTK particle data files into surface meshes using the marching cubes algorithm.

* **`vtk_to_obj_vtk_dir`** *(str, mandatory if `vtk_to_obj=True`)*:
  Path to the directory containing the input VTK files.

* **`vtk_to_obj_out_dir`** *(str, default: same as `vtk_dir`)*:
  Directory where the generated `.obj` mesh files will be saved. If not provided,
  saves the `.obj` files in the same folder as the `.vtk`files

* **`vtk_to_obj_vtk_prefix`** *(str, default: `"PartFluid_"`)*:
  Prefix of the VTK files to be converted. Only files starting with this prefix
  will be processed. This will iterate over `PartFluid_1.vtk`, `PartFluid_2.vtk`
  etc.

* **`vtk_to_obj_particle_radius`** *(float, mandatory if `vtk_to_obj=True`)*:
  Radius of each SPH particle in the input data.

* **`vtk_to_obj_smoothing_length`** *(float, default: `2.0`)*:
  Smoothing length used in the SPH kernel, defined in multiples of the particle radius. The kernel’s compact support radius is twice this value.

* **`vtk_to_obj_cube_size`** *(float, default: `1.0`)*:
  Edge length of the marching cubes grid cells, also in multiples of the particle radius. This defines the resolution of the reconstructed surface.

* **`vtk_to_obj_surface_threshold`** *(float, default: `0.6`)*:
  Iso-surface threshold for the fluid density field, expressed in multiples of the SPH rest density. Determines the surface level used for mesh extraction.

See how you can visualize the results in Blender in the next section [Rendering in Blender](render_in_blender.md).
