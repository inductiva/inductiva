# Getting Started with Inductiva
Welcome to Inductiva, a Python API package that enables running popular [open-source physical system simulation packages]() 
in the Cloud. 

This tutorial provides a step-by-step guide to help you get started with Inductiva 
before running your first simulation project. Should you encounter any errors or 
issues along the way, our [troubleshooting guide]() is readily available to assist 
you in resolving them efficiently.

### Steps We'll Cover:

1. [Install the Inductiva API]()
2. [Request Your API Access Token]()
3. [Verify the Installation with a Test Run]()
4. [Run your First Simulation Example]()
5. [Explore what to read next]()

## Install the Inductiva API with pip
If you've set up the API access token, you're ready to install the 
latest Inductiva package release on PyPI with pip. 

Simply open your terminal and run:

```
pip install --upgrade inductiva
```

Encountering issues? Don’t worry. Head over to our [troubleshooting guide]() to 
work it out.

## Request Your API Access Token
If you don't have a valid API access token, [fill this form](https://docs.google.com/forms/d/e/1FAIpQLSflytIIwzaBE_ZzoRloVm3uTo1OQCH6Cqhw3bhFVnC61s7Wmw/viewform) and request 
your own personal API key.

Once you receive your API token, you can then set the `INDUCTIVA_API_KEY` as an 
environment variable in your terminal:
```
export INDUCTIVA_API_KEY="YOUR_API_KEY"
```

## Verify the Installation with a Test Run

Next, make sure everything is set up correctly with a quick test run.
You have two options to test the installation:
1. through the command line interface (CLI) tool that gets installed
   automatically when the `pip` command runs;
2. programmatically via direct import of the `inductiva` package in a python
   command.

Both options are shown below and both output the same version:

```console
# using the CLI tool ($ is the shell prompt):
$ inductiva --version
inductiva 0.4.2
# programmatically through direct python invocation:
$ python -c 'import inductiva; print(inductiva.__version__)'
0.4.2
```

## Run your First Simulation Example: The Dam Break with Reef3D

In this example, you will use the open-source hydrodynamics REEF3D simulator to
simulate a **2D dam break scenario** where a block of fluid is let to flow under the effect of gravity as follows:

<div align="center">
   <img src="./assets/reef3d-dambreak.gif" alt="REEF3D 2D dambreak simulation" width=500>
</div>

With this first example, you learn that running simulations via Inductiva API is
not much different from running them in your local machine. That's the magic
of it all! 

Running a simulation with API involves preparing a Python script with the following
steps:
1. Preparing all the configuration files for the simulation in a single input folder;
2. Instantiate a simulator object that identifies the simulator you want to use.
In this case, we are going to instantiate the REEF3D simulator object;
3. Launch the simulation with the `run` method of the simulator object and pass
a reference to the input folder. This folder will be uploaded to a user's
remote storage and used in the worker machine that will run the simulation;
4. Wait for the simulation to finish and download the results to your local machine.

To make it simpler, in this example, we will download the input folder with all
the configuration files necessary to run the dam break simulation. These files 
were obtained from the REEF3D tutorials available on their
[GitHub repository](https://github.com/REEF3D/REEF3D/tree/master/Tutorials/REEF3D_CFD/9_1%202D%20Dam%20Break). We
have altered them slightly to reduce the time the simulation takes to run.
Here's how:

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

As the simulation runs information, you will receive information about its progress.
When it finishes, you should see the following folder on your local machine:

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

With the data produced by the simulation on your local machine, you can now post-process
it as you would normally do. For example, you can visualize the contents of the files
in `REEF3D_CFD_VTU` using standard visualization tools such as Paraview.

## What to read next

Learn how to [run one of your own simulation projects]() for the first time through 
the Inductiva API, on the simulator of your choice.

Learn how to [customize the hardware setup]() you use for running your simulations 
through the Inductiva API, and explore the available hardware options to enhance
your project's performance.

If you're looking for inspiration, learn how a group of coastal engineering researchers 
at the University of Porto's Faculty of Engineering have [used the Inductiva API to simulate the most optimal breakwater](https://inductiva.ai/blog/article/scaling-coastal-engineering-projects-inductiva-api) 
to protect some of Portugal's most endangered coastlines.

Troubleshoot installation problems that you might encounter with Inductiva 
API by checking our [troubleshooting guide](#troubleshootingguide), [FAQs](), or
[getting in touch]() with our support team directly.

