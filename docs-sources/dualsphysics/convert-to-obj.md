# Convert VTK Files to Blender-Compatible Formats
To visualize DualSPHysics simulations in Blender, the raw particle data (typically stored in `.vtk` files) 
must first be converted into a mesh format that Blender can interpret, such as `.obj`. This conversion can be easily performed using Inductiva with just a few lines of code.

As an example, consider the [3D dam break](../../multiple_gpus.md) tutorial. In this simulation, we convert the .vtk files to .obj using the following parameters:

To visualize DualSPHysics simulations in Blender, the raw particle data
(typically stored in `.vtk` files) must first be converted into a mesh format
that Blender understands, such as `.obj`. This process can be easily done with
Inductiva with just a few lines of code.

Let's use the [3D Dam Break]() tutorial as an example. In
this simulation we convert the `vtk` files to `obj` with the following parameters:

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

This process uses [Splashsurf](https://github.com/InteractiveComputerGraphics/splashsurf), an open-source Rust-based tool that automatically converts SPH simulation data into high-quality surface meshes saved as `.obj` files.

## Splashsurf for Mesh Generation
Inductiva integrates **Splashsurf** as a post-processing tool, simplifying the conversion process and enabling users to effortlessly import and visualize their simulation results in tools like Blender.

This integration streamlines the workflow from raw simulation data to high-quality visualizations, bridging the gap between numerical data and rendering.

### Post-Processing Options
* **`vtk_to_obj`** *(bool)*:
  Enables the conversion of VTK particle data files into surface meshes using the marching cubes algorithm.

* **`vtk_to_obj_vtk_dir`** *(str, mandatory if `vtk_to_obj=True`)*:
  Path to the directory containing the input VTK files.

* **`vtk_to_obj_out_dir`** *(str, default: same as `vtk_dir`)*:
  Directory where the generated `.obj` mesh files will be saved. Defaults to the same folder as the `.vtk` files if not specified.

* **`vtk_to_obj_vtk_prefix`** *(str, default: `"PartFluid_"`)*:
  Prefix of VTK files to convert. Only files starting with this prefix will be processed (e.g., `PartFluid_1.vtk`, `PartFluid_2.vtk`, etc.)

* **`vtk_to_obj_particle_radius`** *(float, mandatory if `vtk_to_obj=True`)*:
  Radius of each SPH particle in the input data.

* **`vtk_to_obj_smoothing_length`** *(float, default: `2.0`)*:
  The smoothing length used in the SPH kernel, expressed as a multiple of the particle radius. The kernel’s compact support radius is twice this length.

* **`vtk_to_obj_cube_size`** *(float, default: `1.0`)*:
  Edge length of the marching cubes grid cells, also in multiples of the particle radius. Determines the resolution of the reconstructed surface.

* **`vtk_to_obj_surface_threshold`** *(float, default: `0.6`)*:
  Iso-surface threshold for the fluid density field, expressed as a multiple of the SPH rest density. Determines the surface level used for mesh extraction.


Next, learn how to visualize the resulting `.obj` files in Blender in the Rendering in Blender sectio [Rendering in Blender](https://inductiva.ai/guides/dualsphysics/render-in-blender).
