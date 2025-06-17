# The Inductiva-Supercharged Simulation Workflow

Let’s break down each step of running a simulation with Inductiva and see how it transforms your typical workflow for the better!

1️⃣Pick a Cloud Machine -> 2️⃣Pick Your Simulator -> 3️⃣Start Your Simulation -> 4️⃣Stop Your Machine & Analyze Results

## The goal

At the end of this step-by-step guide you will:
💡Master the key commands in the Inductiva simulation workflow.
🐍Have created a python script to run your first simulation.

## Prerequisites

Before starting, ensure the following:

- **Python installed**: Make sure Python (>=3.9) is installed on your system. Check our [System Prep guide](../../systemrequirements/index.md) for more information.
- **Inductiva API Setup**: Ensure Inductiva library is installed. If you haven’t already, go through the installation steps.
For a quick check, see that your API Key is set in [Inductiva's web Console](https://console.inductiva.ai/account/details).

## Create the Python file

Save a new file *my_inductiva_run.py*

- Open your preferred IDE (e.g., VS Code, PyCharm) and save it in your project directory;
- or simply open Notepad++ and save it as .py.
Throughout this guide, you’ll add python commands to this file. Start with the command to import the Inductiva library.

```python
 import inductiva
```

## 1️⃣Pick a Cloud Machine

### List Available Machines

In your computer’s Terminal, use the command

```python
inductiva resources available
```

to display a list of machines, including their types, configurations, and capabilities.
Alternatively, go to the [web Console](https://console.inductiva.ai/machine-groups/instance-types) to check out the list of machines.

Inductiva provides you access to hundreds of cutting-edge cloud machines and empowers you to select the best option for your simulation use-case.
These machines cover a wide range of specifications, and provide you with options at different levels of performance and cost.
Learn how to navigate performance and price options - <a href="pick-cloud-machine.html">Pick a cloud machine for your simulation case</a>

### Allocate your exclusive cloud machine

The machine you selected is the first parameter of the command to create the MachineGroup.

```python
cloud_machine = inductiva.resources.MachineGroup(
 machine_type="c2d-highcpu-112",
 spot=True)
```

💡This instruction will allocate a group of cloud machines exclusively for you - No queue, No waiting time!

🐍 **Python script checkpoint**:
By this time, your *my_inductiva_run.py* file should look like this:

```python
import inductiva

cloud_machine = inductiva.resources.MachineGroup(
 machine_type="c2d-highcpu-112",
 spot=True)
```

## 2️⃣Pick Your Simulator

### List Available Simulators

In your computer’s Terminal, use the command

```python
inductiva simulators ls
```

to list the simulators built in Inductiva API, along with their supported versions.
Alternatively, go to the <a href="https://inductiva.ai/simulators">website</a> to check out the list of simulators.

### Choose a Simulator

Instantiate the selected simulator in your Python script:

```python
reef3d = inductiva.simulators.REEF3D()
```

💡This will initialize the simulator in your MachineGroup.

> 📝 In this guide, we’ll be using reef3d, but here’s a cheat sheet for your preferred simulator:
>
> 1. Go to Inductiva’s simulators webpage.
> 2. Click on table’s cell that matches your simulator.
> 3. In that simulator’s webpage, find the code example and copy the command that initializes that simulator.

🐍 **Python script checkpoint**:
By this time, your *my_inductiva_run.py* file should look like this:

```python
import inductiva

cloud_machine = inductiva.resources.MachineGroup(
 machine_type="c2d-highcpu-112",
 spot=True)

reef3d = inductiva.simulators.REEF3D()
```

## 3️⃣Start Your Simulation

Once your simulator is set on your MachineGroup, start your simulation and wait for it to finish.

```python
task = reef3d.run(
  input_dir="/path/to/my/input/files",
  on=cloud_machine)
 task.wait()
```

💡This will:

- Run the simulator,
- With the input files stored in the indicated directory path,
- On your MachineGroup.

## 4️⃣Stop Your Machine & Get the Results

When your simulation is finished, terminate your cloud machine and download your simulation results.

```python
 cloud_machine.terminate()
 task.download_outputs()
```

💡An idle machine still costs you credits. Terminate it after the task is finished to save you compute costs.
💡If you don’t do it explicitly, Inductiva will do it for you after a while.

🐍 **Python script checkpoint**:
By this time, your *my_inductiva_run.py* file should look like this:

```python
import inductiva

cloud_machine = inductiva.resources.MachineGroup(
 machine_type="c2d-highcpu-112",
 spot=True)

reef3d = inductiva.simulators.REEF3D()

task = reef3d.run(
  input_dir="/path/to/my/input/files",
  on=cloud_machine)
 task.wait()

cloud_machine.terminate()
task.download_outputs()
```

**Your Python script is ready! Let’s run your simulation**
Run your *my_inductiva_run.py* file, using your preferred option:
On your IDE: run the project.
On the computer’s Terminal:
   navigate to the location of the .py file:

```python
cd path/to/your/.py/file
```

   and run the file:

```python
python my_inductiva_run.py
```
