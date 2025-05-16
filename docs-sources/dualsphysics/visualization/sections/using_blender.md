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

<p align="center"><img src="./_static/dam_break_elastic.gif" alt="Visualization created with Blender." width="700"></p>


```{toctree}
:hidden:
Convert vtk to obj <convert_to_obj.md>
Using Blender <render_in_blender.md>
```