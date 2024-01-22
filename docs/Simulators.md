# Simulators

Inductiva API has available several open-source simulators ready to use. Users who are familiar with the simulators can easily start running simulations with their previously prepared simulation configuration files. In this way, they can take advantage of performant hardware to speed up their simulation and exploration.
No installation or management of the simulators is required, no need to worry about hardware or software dependencies. Inductiva API takes care of all of that for you.

The simulators currently available are all open-source and have their dedicated documentation:
- [SPlisHSPlasH](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH)
- [DualSPHysics](https://dual.sphysics.org/)
- [OpenFOAM](https://www.openfoam.com/)
- [SWASH](https://swash.sourceforge.io/)
- [XBeach](https://oss.deltares.nl/web/xbeach/)
- [Reef3D](https://reef3d.wordpress.com/)
- [GROMACS](https://www.gromacs.org/)
- [FDS](https://pages.nist.gov/fds-smv/)

Check the documentation of each simulator to learn more on how to configure them. Here, we highlight how to use your preferred simulator via Inductiva API and how to scale your simulations with simplicity. There is a general structure that all simulators follow, and then there are some specificities for each simulator.

## Simulators via Inductiva API

To run a simulation, prepare the configuration files for your desired simulator and place them in a designated folder. The simulator will use this folder to perform the simulation.

### Example

For a first simulation, we run a coastal dynamics simulation with SWASH.

```python
import inductiva

input_dir = inductiva.utils.files.download_from_url(
    "https://storage.googleapis.com/inductiva-api-demo-files/"
    "swash-input-example.zip")

# Initialize the Simulator
swash = inductiva.simulators.SWASH()

# Run simulation with config files in the input directory
task = swash.run(input_dir=input_dir, 
                 sim_config_filename="input.sws")

task.download_outputs()
```

And that's it! Your simulation is now running in the cloud, and you have a `task` object that allows you to manage it. You can check its status with `task.get_status()`, wait for it to finish with `task.wait()`, and download the results with `task.download_outputs()`.

With Inductiva API you don't have immediate access to visualization tools of these simulators. However, you can download the results and use the visualization tools of your choice. 

## Large-scale simulations

In all of the examples above, your simulations will run on our default set of machines, which are available for all users to use. These machines are not the most performant and are mostly useful for demo and testing purposes. To run large-scale simulations you have the option to choose your own dedicated machine group, where only your simulations will run.

For instance, let's use DualSPHysics as an example again. Let's launch a machine with 15 CPU physical cores and 120 GB of RAM and run our simulation there. Learn further on how to manage the computational resources of Inductiva API [here](https://github.com/inductiva/inductiva/blob/main/inductiva/resources/README.md).

### Example

```python
import inductiva

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