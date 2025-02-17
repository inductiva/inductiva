# How to access simulation output while it is running

This document provides a step-by-step guide on how to monitor and access simulation output while the simulation is still running. By following this guide, you will learn how to start your simulation, list the output files, and tail specific files for real-time updates. These tools are useful for debugging and tracking the progress of your simulation.


## Start your simulation
First, you need to start a `Task` using the [Inductiva API](https://inductiva.ai/). If you are not familiar with the API we suggest you follow [this tutorial](https://docs.inductiva.ai/en/latest/intro_to_api/tasks.html) to learn more about launching simulations.
To show the process of accessing simulation output in real-time, we will run the example [DualSPhysics](https://tutorials.inductiva.ai/simulators/DualSPHysics.html) simulation:

```python
import inductiva

# Instantiate machine group
cloud_machine = inductiva.resources.MachineGroup(
    provider="GCP",
    machine_type="c2-standard-4")

# Download the configuration files into a folder
input_dir = inductiva.utils.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "dualsphysics-input-example.zip",
    unzip=True)

# Initialize the Simulator
dualsphysics = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = dualsphysics.run(
    input_dir=input_dir,
    shell_script="run.sh",
    on=cloud_machine)

task.wait()
cloud_machine.terminate()

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
config_out/
├── config_out/CfgInit_Normals.vtk
├── config_out/CfgInOut_DomainBox.vtk
├── config_out/config__Dp.vtk
├── config_out/config.xml
├── config_out/config_Bound.vtk
├── config_out/CfgInit_NormalsGhost.vtk
├── config_out/Run.out
├── config_out/data/
│   ├── config_out/data/Part_0004.bi4
│   ├── config_out/data/Part_0000.bi4
│   ├── config_out/data/Part_0002.bi4
│   ├── config_out/data/Part_0003.bi4
│   ├── config_out/data/Part_0005.bi4
│   ├── config_out/data/Part_Head.ibi4
│   ├── config_out/data/PartOut_000.obi4
│   ├── config_out/data/PartInfo.ibi4
│   └── config_out/data/Part_0001.bi4
├── config_out/config_Fluid.vtk
├── config_out/config_All.vtk
├── config_out/CfgInit_MapCells.vtk
├── config_out/config.bi4
├── config_out/CfgInit_Domain.vtk
├── config_out/config.out
├── config_out/CfgInOut_DomainReal.vtk
└── config_out/config_MkCells.vtk
config.xml
stderr.txt
run.sh
config_Def.xml
stdout.txt
```

## Tail specific file
In order to read the last lines of an output file, we introduce the tail command, `inductiva task tail <task_id> <filename>`. This command allows you to track updates as they are written to a file: Make sure that the file you are trying to tail is a text file. 

For example 
```bash
$ inductiva tasks tail jjqwbn2yi5kwk5rv6qs2nth21 config_out/Run.out
```
Prints out the last 10 lines of the `config_out/Run.out` file:

```bash
┌ (last 10 lines from config_out/Run.out)
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
$ inductiva tasks tail jjqwbn2yi5kwk5rv6qs2nth21 config_out/Run.out --lines 20
```

Prints out the last 20 lines of the file:

```bash
┌ (last 20 lines from config_out/Run.out)

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
$ inductiva tasks tail jjqwbn2yi5kwk5rv6qs2nth21 config_out/Run.out -f
```
To stop monitoring, press `Ctrl+C`

This command is especially useful for tracking the progress of the simulation, debugging issues, and identifying potential problems early, allowing you to terminate your task preemptively if needed.