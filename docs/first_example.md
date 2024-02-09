# Run your first simulation

In this example, you will use the open-source [hydrodynamics REEF3D simulator](https://github.com/REEF3D/REEF3D) to
simulate a **2D dam break scenario**. This involves a block of fluid released to 
flow under the influence of gravity, as shown in the video below. As the simulation
progresses, you'll find that running simulations through the Inductiva API is not 
much different from running them in your local machine.

<div align="center">
   <img src="./_static/reef3d-dambreak-fullscreen.gif" alt="REEF3D 2D dambreak simulation">
</div>


To run a simulation via the API, you first need to prepare a Python script with 
the following steps:

1. **Prepare and gather all the necessary configuration files for the simulation 
into one input folder**. To simplify our example, we will provide you with an input folder 
containing all the necessary configuration files for the dam break simulation, 
slightly modified from the original REEF3D tutorials found on their [GitHub repository](https://github.com/REEF3D/REEF3D/tree/master/Tutorials/REEF3D_CFD/9_1%202D%20Dam%20Break) to reduce
the simulation run time;

2. **Instantiate a simulator object that identifies the simulator you want to use**, 
here exemplified by instantiating a REEF3D simulator object;

3. **Launch the simulation using the `run` method of the simulator object, passing
a reference to the input folder**. This folder gets uploaded to your remote storage 
and accessed by the worker machine executing the simulation;

4. After the simulation completes, **download the results back to your local machine.**

Following the steps outlined above, here is the Python script you'll use to run 
the dam break simulation via the Inductiva API. Give it a go:

```python
import inductiva

# 1 - Download the configuration files for a REEF3D simulation. This folder
# contains the control.txt and ctrl.txt files that configure the parameters of
# the mesh and the simulation, respectively.
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "reef3d-dambreak-example.zip", unzip=True)

# 2 - Initialize the REEF3D simulator object
simulator = inductiva.simulators.REEF3D()

# 3 - Launch the simulation with the downloaded input folder. This will return
# a task object that can be used to monitor the task and download the outputs.
task = simulator.run(input_dir=input_dir)

# 4 - Wait for the simulation to finish and download the outputs to a default
# folder in your local machine named `example_simulation`.
task.wait()
task.download_outputs(output_dir="example_simulation")
```
As the simulation progresses, you'll receive updates on its status. When it finishes, 
you should find the following folder on your local machine:

```bash
inductiva_output/example_simulation
|
|- DIVEMesh_Decomp
|- DIVEMesh_Log
|- DIVEMesh_Paraview 
|- REEF3D_CFD_VTU
|- REEF3D_Log
|- control.txt
|- ctrl.txt
|- grid-000001.dat
|- grid-000002.dat
|- stderr.txt
|- stdout.txt
```
Once the simulation data is on your local machine, you're ready to proceed with 
post-processing just as you normally would. For example, you can visualize the 
contents of the files in `REEF3D_CFD_VTU` using open-source visualization tools 
such as [Paraview](https://www.paraview.org/download/).

## What to read next

Learn about the [shared and dedicated resources](./introduction/computational_resources_overview.md) 
you can use for running your simulations through the Inductiva API, and explore 
the available hardware options to enhance
your project's performance.

If you're looking for inspiration, learn how a group of coastal engineering researchers 
at the University of Porto's Faculty of Engineering have [used the Inductiva API to simulate the most optimal breakwater](https://inductiva.ai/blog/article/scaling-coastal-engineering-projects-inductiva-api) 
to protect some of Portugal's most endangered coastlines.