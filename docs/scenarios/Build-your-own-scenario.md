# Build your own scenario

**Inductiva API** allows users to build their own scenarios and run simulations on them. In this way, users will be able to iterate over their parameters of choice as they wish, without needing to create input files for each of them. The main steps to build a scenario are:

- Template your input files with your parameters of choice.
- Create your scenario class that inherits from `inductiva.scenarios.Scenario`.
- Run your simulations.

## Template your input files

The first step to building your scenario is to template your input files. In this way, users can configure certain parameters to be more easily changeable through Python code. This is done by substituting the parameters of interest with a template variable like `{{ template_arg }}`. 

Let's follow the steps to template one input file of the SplishSplash simulator. In this case, we want to template the following input file  `input.xml` for SplishSplash:

<p align="center">
  <img src="assets/splishsplash_input_file.png" alt="SplishSplash input file" width="550" height="450">
</p>

In this example, we will template the first three arguments presented: the simulation time, the time step size and the particle radius. To do so, we substitute the values presented in the above file with the following template variables: `{{ simulation_time }}`, `{{ time_step }}` and `{{ particle_radius }}`. The templated file looks as follows:

<p align="center">
  <img src="assets/template.png" alt="Template file" width="550" height="450">
</p>

All files that are templated just have the append to their name in the end `.jinja`. In this case, the file is named `input.xml.jinja`.
That is it! Now the rest follows in Python code! 

## Create your scenario class

With your template arguments in mind, you can create your scenario starting from the base scenario `inductiva.scenarios.Scenario`. Based on the template above, we create a simple scenario named `MyScenario` to explain the process. 

### Prepare template directory
The first step is to set up the template directory where the template files are stored. This directory can be placed and named as you wish. However, it must have the following structure:

```
|- template_files_dir
|    |- sim_config_files
|    |    |- input_dir_1
|    |    |    |- input_file_1.jinja
|    |    |- input_file_2.jinja
|    |- commands.json
```

Following this structure, we will be able to process everything correctly. The `commands.json` file is a file that contains the commands to run the simulation. However, only simulators with multiple commands use it, like OpenFOAM. This file can also be templated in the same manner as the other files.

So for our example, we save our template directory in the same directory as the file running our scenario and it will look as follows:
```
|- myscenario_template_dir
|    |- sim_config_files
|    |    |- input.xml.jinja
```

### Create the scenario class

With the template and arguments ready, we are prepared to start writing the scenario class. 

Let's start by importing the necessary modules and setting up the valid simulators and the template directory. 

```python
from inductiva import scenarios, tasks, simulators, resources

class MyScenario(scenarios.Scenario):
    """MyScenario class."""

    valid_simulators = [simulators.SplishSplash]
    template_files_dir = "myscenario_template_dir"
```

The `valid_simulators` attribute is a list of simulators that can be used in this scenario. The `template_files_dir` attribute is the path to the template directory we created above.

The next step is to initialise the scenario. Here, you can set the template parameters you wish - either all of them or only some of them. For structural reasons, we tend to set the parameters more concerned with the nature of the scenario. In this case, we pass the particle radius.

```python
    def __init__(self, particle_radius: float = 0.02):
        """Initialize a MyScenario object.
        
        Args:
            particle_radius: Radius of the fluid particles, in meters.
        """

        self.params["particle_radius"] = particle_radius
```

Then, to finish you set the simulate method. Usually, we decide to pass here the template parameters 
more concerned with the simulation configuration. In this case, we pass the simulation time and the time step size. 

**Remark**: There are some input arguments which are default and should always be added. These are the `simulator`, `machine_group` and `storage_dir` arguments. The `simulator` argument allows the configuration of the simulator to be used in the simulation. The `machine_group` argument is the machine group to be used in the simulation. The `storage_dir` argument is the remote directory path where the simulation results will be stored.

```python
    def simulate(
        self,
        simulator: simulators.Simulator = simulators.SplishSplash(),
        machine_group: Optional[resources.MachineGroup] = None,
        storage_dir: Optional[str] = "",
        time_step: float = 0.01,
        simulation_time: float = 1):
        """Simulates the scenario."""

        # Set the template args
        self.params["time_step"] = time_step
        self.params["simulation_time"] = simulation_time

        # Create the task
        # The first three arguments: simulator, machine_group and
        # storage_dir are obligatory to be passed to the base class.
        # The `sim_config_filename` argument is the name of the input
        # file and is required for SplishSplash.
        task = super().simulate(
            simulator=simulator,
            machine_group=machine_group,
            storage_dir=storage_dir,
            sim_config_filename="input.xml"
        )

        return task
```

And that is mostly it! In this way, you can use your scenario to iterate over the parameters you templated more easily. 

There are a few extra ways to set up your scenario, and we dig into them below. Feel free to notify us of any questions or suggestions you may have!

### Extra functionalities

Currently, there are two extra functions users can define in their class to configure further the scenario. These are the `config_params` and `add_extra_input_files` functions.

The `config_params` allows you to further configure the parameters of the simulation. For example, you can set a parameter in the initialisation of the scenario and then change it in the `config_params` function. 

A simple example would be to set the resolution argument in the `__init__`. In the template file, the argument that controls the resolution is the particle radius, therefore, we can use the `config_params` to set that. The purpose of this function is to not clutter the initialisation and simulate methods, however, what is done here, could be done there.

The `add_extra_input_files` permit users to add extra files as inputs to the simulation. The main goal is to allow input files to also be parameterizable. For example, adding objects to the simulation. In this way, users should set these file in the appropriate place for the simulation to run with this function. Here, the file should be placed in the `input_dir` which will only contain the simulation files.
