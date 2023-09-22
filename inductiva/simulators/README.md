# Simulators

Inductiva API has available several open-source simulators ready to use. Users that are familiar with the simulators can easily start running simulations with their previously prepared simulation configuration files. In this way, they can take advantage of performant hardware to speed up their simulation and exploration.
No installation or management of the simulators is required, no need to worry about hardware or software dependencies. Inductiva API takes care of all of that for you.

The simulators currently available are all open-source and have their own dedicated documentation:
- [SPlisHSPlasH](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH)
- [DualSPHysics](https://github.com/DualSPHysics/DualSPHysics)
- [OpenFOAM](https://www.openfoam.com/)
- [SWASH](https://swash.sourceforge.io/)
- [XBeach](https://oss.deltares.nl/web/xbeach/)
- [GROMACS](https://www.gromacs.org/)

Check the documentation of each simulator to learn more on how to configure them. Here, we highlight how to use your preferred simulator via Inductiva API and how to scale your simulations with simplicity. There is a general structure that all simulators follow, and then there are some specificities for each simulator.

## Simulators via Inductiva API

To run a simulation, prepare the configuration files for your desired simulator and place them in a designated folder. The simulator will use this folder to perform the simulation.

Here, we present an example using DualSPHysics - a Smoothed-Particle Hydrodynamics simulator - and show how it looks in practice. To simplify we provide the configuration files for this example.

### Example

This follows a classical example in CFD with flow around a cylinder.

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Download the configuration files into a folder
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "dualsph-flow-cylinder.zip"
)

# Initialize the Simulator
simulator = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = simulator.run(input_dir=input_dir)
```

And that's it! Your simulation is now running in the cloud, and you have a `task` object that allows you to manage it. You can check its status with `task.get_status()`, wait for it to finish with `task.wait()`, and download the results with `task.get_output()`.

With Inductiva API you don't have immediate access to visualization tools of these simulators. However, you can download the results and use the visualization tools of your choice. For example, for DualSPHysics you can use [Paraview](https://www.paraview.org/) to visualize the results.

## DualSPHysics simulator

DualSPHysics is a Smoothed-Particle Hydrodynamics (SPH) simulator. The simulator is usually configured by a single file with the extension `.xml`. This file contains all the information about the simulation, including the geometry, the physical properties of the fluids, the boundary conditions, the numerical parameters, and the output files. Sometimes the configuration can also use extra geometry files. 

The default name used for the configuration file is `config.xml`. However, you can use any name for the configuration file when running the simulation, e.g., `my_simulation.xml`:

```python
task = simulator.run(input_dir=input_dir, sim_config_file="my_simulation.xml")
```

## SPlisHSPlasH simulator

SPlisHSPlasH is a Smoothed-Particle Hydrodynamics (SPH) simulator that covers a wide-range of applications. The simulator is usually configured by a single file with the extension `.json`. This file contains all the information about the simulation, including the geometry, the physical properties of the fluids, the boundary conditions, the numerical parameters, and the output files. Sometimes the configuration can also use extra geometry files.

To run a simulation you can follow exactly the same structure as above for DualSPHysics.

### Example

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Set simulation input directory
input_dir = "splishsplash_example"

# Initialize the Simulator
simulator = inductiva.simulators.SplishSplash()

# Run simulation with config files in the input directory
task = simulator.run(input_dir=input_dir,
                     sim_config_file="config.json")
```

## OpenFOAM simulator

OpenFOAM is a Finite Volume method for CFD simulations with a wide-range of applications across several areas of engineering and science. OpenFOAM has an extensive
range of features to solve anything from complex fluid flows involving chemical reactions, turbulence and heat transfer, to solid dynamics and electromagnetics.

A single simulation via Inductiva API comprises several steps done via OpenFOAM - e.g., partitioning the domain, meshing, solvers and post-processing. Hence, to configure a simulation for OpenFOAM the user will need a set of configuration files that are organizaed in three folders:
- `time`: containing individual files of data for particular fields, like initial values and boundary conditions that the user must specify. For example, for an initial condition at $t=0$ the initial conditions will be stored in the directory `0`.
- `constant`: contains files that describe the objects in the simulation and the physical properties for the application we are concerned.
- `system`: contains all of the files that describe the simulation, including the solvers, the numerical parameters, and the output files. It must contain at least 3 files: `controlDict` where run control parameters are set including start/end time, time step and parameters for data output; `fvSchemes` where discretisation schemes used in the solution may be selected at run-time; and `fvSolution` where the equation solvers, tolerances and other algorithm controls are set for the run.

All of these folders should be inside an input directory. Finally, to run a simulation the user needs to configure a list of dictionaries specifying the commands they want to execute on the backend. Below, we run the [motorbike tutorial](https://github.com/OpenFOAM/OpenFOAM-8/tree/master/tutorials/incompressible/simpleFoam/motorBike) from OpenFOAM and show how this is done in practice.

### Example 

````python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Set simulation input directory
input_dir = "motorbike_tutorial"

# Set the simulation commands
[{"cmd": "surfaceFeatures", "prompts": []}, {"cmd": "blockMesh", "prompts": []},
{"cmd": "snappyHexMesh", "prompts": []}, {"cmd": "potentialFoam", "prompts": []},
{"cmd": "simpleFoam", "prompts": []}]

# Initialize the Simulator
simulator = inductiva.simulators.OpenFOAM()

# Run simulation with config files in the input directory
task = simulator.run(input_dir=input_dir, commands=commands)
````

## SWASH simulator

SWASH is a simulator that solves the shallow water equations and is used to simulate waves and currents in coastal waters and harbours, long waves in coastal regions and tidal inlets, and rapidly-varied flows around coastal structures. The simulator is configured using a single file with the `.sws` extension, and additional files containing information about the domain and the ocean floor, such as a bathymetry file with a `.bot` extension, are necessary for the simulation to run.

### Example

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Set simulation input directory
input_dir = "swash_example"

# Initialize the Simulator
simulator = inductiva.simulators.SWASH()

# Run simulation with config files in the input directory
task = simulator.run(input_dir=input_dir, sim_config_filename="input.sws")
```

## XBeach simulator

XBeach is a simulator with a two-dimensional model for wave propagation, sediment transport and morphological changes of the nearshore area. The simulator is configured with a `params.txt` file that contains grid and bathymetry info, wave input, flow input, morphological input, etc. in the form of keyword/value pairs. If a `params.txt` cannot be found then XBeach will not run. Other files are used to configure the grid and bathymetry profile, like `bed.dep` for example, other files with extra information that can be used inside the `params.txt` to configure the simulator further.

### Example

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Set simulation input directory
input_dir = "swash_example"

# Initialize the Simulator
simulator = inductiva.simulators.XBeach()

# Run simulation with config files in the input directory
task = simulator.run(input_dir=input_dir)
```

## GROMACS simulator

GROMACS is a versatile simulator to perform molecular dynamics simulations. It is primarily designed for biochemical molecules like proteins, lipids and nucleic acids that have a lot of complicated bonded interactions, but since GROMACS is extremely fast at calculating the nonbonded interactions (that usually dominate simulations) many groups are also using it for research on non-biological systems, e.g. polymers and fluid dynamics.

A single simulation of GROMACS via Inductiva API can comprise several steps - e.g., preparing the molecules, minimizing the energy of the system, running the simulation and post-processing. Hence, to configure a simulation of GROMACS the user may require several files. Moreover, GROMACS has specific commands to run certain tasks that already use the files in your input directory. 

Below, we run the [protein solvation scenario]() with GROMACS and show how this is done in practice.

### Example

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Set simulation input directory
input_dir = "molecules_water_box"

# Set the simulation commands
commands = [
    {"cmd": "gmx solvate -cs tip4p -box {{ box_size }} -o conf.gro -p topol.top", "prompts": []},
    {"cmd": "gmx grompp -f energy_minimization.mdp -o min.tpr -pp min.top -po min.mdp -c conf.gro -p topol.top", "prompts": []},
    {"cmd": "gmx mdrun -s min.tpr -o min.trr -c min.gro -e min.edr -g min.log", "prompts": []},
    {"cmd": "gmx grompp -f positions_decorrelation.mdp -o decorr.tpr -pp decorr.top -po decorr.mdp -c min.gro", "prompts": []},
    {"cmd": "gmx mdrun -s decorr.tpr -o decorr.trr -x  -c decorr.gro -e decorr.edr -g decorr.log", "prompts": []},
    {"cmd": "gmx grompp -f simulation.mdp -o eql.tpr -pp eql.top -po eql.mdp -c decorr.gro", "prompts": []},
    {"cmd": "gmx mdrun -s eql.tpr -o eql.trr -x eql.xtc -c eql.gro -e eql.edr -g eql.log", "prompts": []}
]

# Initialize the Simulator
simulator = inductiva.simulators.GROMACS()

# Run simulation with config files in the input directory
task = simulator.run(input_dir=input_dir, commands=commands)
```

## Large-scale simulations

In all of the examples above, your simulations will run on our default set of machines, which are available for all users to use. These machines are not the most performant and are mostly useful for demo and testing purposes. To run large-scale simulations you have the option to choose your own dedicated machine group, where only your simulations will run.

For instance, let's use DualSPHysics as an example again. Let's launch a machine with 15 CPU physical cores and 120 GB of RAM and run our simulation there. Learn further on how to manage the computational resources of Inductiva API [here](https://github.com/inductiva/inductiva/blob/main/inductiva/resources/README.md).

### Example

```python
import inductiva

inductiva.api_key = "YOUR_API_KEY"

# Initialize your machines
my_machine_group = inductiva.resources.MachineGroup(
    machine_type="c2d-standard-30,
    num_machines=1,
    disk_size_gb=50
)

# Start the machines
my_machine_group.start()

# Download the configuration files into a folder
input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "dualsph-flow-cylinder.zip"
)

# Initialize the Simulator
simulator = inductiva.simulators.DualSPHysics()

# Run simulation with config files in the input directory
task = simulator.run(input_dir=input_dir, machine_group=my_machine_group)

# Wait for the simulation to finish
task.wait()

# Terminate the machine
my_machine_group.terminate()
```