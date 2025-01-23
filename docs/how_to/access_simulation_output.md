# How to access simulation output while it is running

This document provides a step-by-step guide on how to monitor and access simulation output while the simulation is still running. By following this guide, you will learn how to start your simulation, list the output files, and tail specific files for real-time updates. These tools are useful for debugging and tracking the progress of your simulation.


## Start your simulation
First, you need to start a `Task` using the [Inductiva API](https://inductiva.ai/). If you are not familiar with the API we suggest you follow [this tutorial](https://docs.inductiva.ai/en/latest/intro_to_api/tasks.html) to learn more about launching simulations.
To show the , we will run the example [DualSPhysics](https://tutorials.inductiva.ai/simulators/DualSPHysics.html) simulation:

```python
import inductiva

# Instantiate machine group
machine_group = inductiva.resources.MachineGroup("c2-standard-4")
machine_group.start()

# Download the configuration files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "dualsphysics-input-example.zip",
    unzip=True)

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = dualsphysics.run(input_dir=input_dir,
                        shell_script="run.sh",
                        on=machine_group)

task.wait()
machine_group.terminate()

task.download_outputs()

task.print_summary()
```

## List Files
We introduce the command inductiva task list-files <task_id>, which lists all files recursively in the output directory of the specified task.

For example:

```bash
$ inductiva tasks list-files jjqwbn2yi5kwk5rv6qs2nth21
```

Produces an output similar to:

```bash
Directory contents:
[FILE] stderr.txt
[FILE] flow_cylinder.xml
[FILE] flow_cylinder__Dp.vtk
[FILE] stdout.txt
[FILE] flow_cylinder_Fluid.vtk
[FILE] flow_cylinder.bi4
[DIR] flow_cylinder
  [FILE] Run.out
  [FILE] CfgInit_NormalsGhost.vtk
  [FILE] CfgInit_Domain.vtk
  [FILE] CfgInit_MapCells.vtk
  [FILE] CfgInOut_DomainBox.vtk
  [DIR] data
    [FILE] Part_0004.bi4
    [FILE] Part_0000.bi4
    [FILE] PartOut_000.obi4
    [FILE] Part_0001.bi4
    [FILE] Part_Head.ibi4
    [FILE] Part_0003.bi4
    [FILE] Part_0005.bi4
    [FILE] PartInfo.ibi4
    [FILE] Part_0002.bi4
  [FILE] CfgInit_Normals.vtk
  [FILE] CfgInOut_DomainReal.vtk
[FILE] flow_cylinder_Bound.vtk
[FILE] config.xml
[FILE] flow_cylinder.out
[FILE] flow_cylinder_All.vtk
[FILE] flow_cylinder_MkCells.vtk
```

## Tail specific file
In order to read the last lines of an output file, we introduce the tail command, `inductiva task tail <task_id> <filename>`. This command allows you to track updates as they are written to a file: Make sure that the file you are trying to tail is a text file. 

For example 
```bash
$ inductiva tasks tail jjqwbn2yi5kwk5rv6qs2nth21 flow_cylinder/Run.out
```
Prints out the last 10 lines of the `flow_cylinder/Run.outflow_cylinder/Run.out` file:

```bash
┌ (last 10 lines from flow_cylinder/Run.out)

│Part_0024      0.600229          1535       63      87.02  16-12-2024 16:48:49
│  Particles new: 99 (total new: 3762)  -  Current np: 13599
│  Particles out: 9  (total out: 264)
│Part_0025      0.625266          1598       63      86.15  16-12-2024 16:48:50
│  Particles new: 99 (total new: 3861)  -  Current np: 13506
│  Particles out: 6  (total out: 270)
│Part_0026      0.650338          1660       62      87.60  16-12-2024 16:48:50
│  Particles new: 198 (total new: 4059)  -  Current np: 13601
│  Particles out: 6  (total out: 276)
│
└

```

If you wish to configure the number of lines of the output you can use the argument `--lines`.

The command:
```bash
$ inductiva tasks tail jjqwbn2yi5kwk5rv6qs2nth21 flow_cylinder/Run.out --lines 20
```

Prints out the last 20 lines of the file:

```bash
┌ (last 20 lines from flow_cylinder/Run.out)

│  Particles out: 11  (total out: 577)
│Part_0073      1.825354          4537       60      27.15  23-12-2024 10:47:12
│  Particles new: 99 (total new: 9801)  -  Current np: 13338
│  Particles out: 7  (total out: 584)
│Part_0074      1.850040          4597       60      24.95  23-12-2024 10:47:12
│  Particles new: 198 (total new: 9999)  -  Current np: 13418
│  Particles out: 4  (total out: 588)
│Part_0075      1.875261          4659       62      23.77  23-12-2024 10:47:12
│  Particles new: 99 (total new: 10098)  -  Current np: 13392
│  Particles out: 4  (total out: 592)
│Part_0076      1.900121          4720       61      24.02  23-12-2024 10:47:11
│  Particles new: 99 (total new: 10197)  -  Current np: 13370
│  Particles out: 6  (total out: 598)
│Part_0077      1.925168          4781       61      24.95  23-12-2024 10:47:12
│  Particles new: 99 (total new: 10296)  -  Current np: 13338
│  Particles out: 4  (total out: 602)
│Part_0078      1.950401          4843       62      35.30  23-12-2024 10:47:12
│  Particles new: 198 (total new: 10494)  -  Current np: 13414
│  Particles out: 1  (total out: 603)
│
└
```

To monitor file updates in real-time, use the -f argument. Run the following command:
```bash
$ inductiva tasks tail jjqwbn2yi5kwk5rv6qs2nth21 flow_cylinder/Run.out -f
```
To stop monitoring, press `Ctrl+C`

This command is especially useful for tracking the progress of the simulation, debugging issues, and identifying potential problems early, allowing you to terminate your task preemptively if needed.